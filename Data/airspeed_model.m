% Derive a counts -> airspeed calibration for the vane sensor against the
% VDAS pitot, and characterise encoder noise from the zero-airspeed period.
%
% Physics: drag on the vane ~ q = 0.5*rho*v^2, balanced by a restoring
% torque ~ theta, so theta ~ q and v ~ sqrt(theta). The fit is done in the
% dynamic-pressure domain (q vs theta), where the pitot's 1 Pa integer
% quantisation is uniform and the model is linear in its parameters.
%
% Board pipeline modelled: 2 kHz sampling -> 10-sample block mean -> 200 Hz
% serial output -> moving average of up to 10 samples (50 ms).

clear; close all;

%% Config
rootDir  = fileparts(mfilename('fullpath'));
logDir   = fullfile(rootDir, 'logs');
figDir   = fullfile(rootDir, 'figures');
runs     = {'run_20261029_125514', 'run_20261029_125642'};
labels   = {'Ramp', 'Step'};
ENC_BITS = 14;
ENC_RES  = 2^ENC_BITS;              % counts per revolution
DEG_PER_COUNT = 360 / ENC_RES;
WIN_S    = 0.25;                    % averaging window matched to pitot rate
Q_MIN_FIT = 3;                      % Pa; below this the 1 Pa pitot quantisation dominates
DECIM    = 10;                      % on-board block average, 2 kHz -> 200 Hz
OUT_AVG  = [1 2 5 10];              % candidate moving-average lengths at 200 Hz
% Steady-flow segments for the spectral analysis: [run, t0, t1] in s from
% run start, picked from the pitot trace (step run plateaus + ramp top)
PLATEAUS = [2 23 28; 2 34 63; 2 68.5 75; 2 77 84; 2 86 93; 1 47 56];

if ~exist(figDir, 'dir'), mkdir(figDir); end
set(groot, 'defaultFigureVisible', 'off');
cleanupVis = onCleanup(@() set(groot, 'defaultFigureVisible', 'on'));

%% Load runs and rebuild the sensor time base
% The logged 2 kHz timestamps are host-assigned in batches (they step back
% by up to ~15 ms on resync), so time is reconstructed from the sequence
% counter with a linear fit against the host clock.
R = struct([]);
for k = 1:numel(runs)
    s = readtable(fullfile(logDir, [runs{k} '_counts_student.csv']));
    g = readtable(fullfile(logDir, [runs{k} '_counts_vdas.csv']));

    n  = double(s.seq - s.seq(1));
    p  = polyfit(n, s.timestamp_epoch - s.timestamp_epoch(1), 1);
    ts = s.timestamp_epoch(1) + polyval(p, n);
    jit = s.timestamp_epoch - ts;

    % Unwrap the 14-bit count in case the vane crosses the encoder rollover
    c = double(s.count);
    c = c - ENC_RES * cumsum([0; round(diff(c) / ENC_RES)]);

    % Zero period: before the 1 s moving mean leaves the initial rest value
    cInit = median(c(1:1000));
    iz = find(abs(movmean(c, 2000) - cInit) > 8, 1) - 2000;
    c0 = mean(c(1:iz));

    % Pitot: pa_raw is integer Pa with negative = positive dynamic pressure
    q   = -double(g.pa_raw);
    pos = g.pa_raw ~= 0;
    rho = median(2 * abs(g.pa_raw(pos)) ./ g.value_raw_mps(pos).^2);

    R(k).name = runs{k};   R(k).label = labels{k};
    R(k).ts = ts;          R(k).c = c;   R(k).fs = 1 / p(1);
    R(k).jitter = jit;     R(k).iz = iz; R(k).c0 = c0;
    R(k).cEnd = median(c(end-2000:end));
    R(k).tg = g.timestamp_epoch;  R(k).q = q;  R(k).vgt = g.value_raw_mps .* sign(q);
    R(k).rho = rho;
end
rho = mean([R.rho]);

fprintf('\n===== TIME BASE =====\n');
for k = 1:numel(R)
    fprintf('%-5s: effective sample rate %.2f Hz (nominal 2000), host timestamp jitter std %.2f ms, max %.1f ms\n', ...
        R(k).label, R(k).fs, 1e3*std(R(k).jitter), 1e3*max(abs(R(k).jitter)));
end
fprintf('Air density implied by VDAS v = sqrt(2q/rho): rho = %.4f kg/m^3\n', rho);

%% Encoder noise at zero airspeed
fprintf('\n===== ENCODER NOISE (zero airspeed) =====\n');
fprintf('Resolution: %d-bit, 1 count = %.4f deg; ideal quantisation noise = 1/sqrt(12) = %.3f counts\n', ...
    ENC_BITS, DEG_PER_COUNT, 1/sqrt(12));
noise = struct([]);
for k = 1:numel(R)
    z  = R(k).c(1:R(k).iz);
    zd = z - mean(z);
    lin = polyfit((0:numel(z)-1)'/R(k).fs, z, 1);
    [pxx, f] = pwelch(zd, hann(4096), 2048, 4096, R(k).fs);
    [tau, adev] = allan_dev(zd, R(k).fs);
    ac = xcorr(zd, 20, 'coeff'); ac = ac(21:end);

    noise(k).std = std(z);  noise(k).p2p = range(z);  noise(k).dur = numel(z) / R(k).fs;
    noise(k).drift = lin(1);  noise(k).f = f;  noise(k).pxx = pxx;
    noise(k).tau = tau;  noise(k).adev = adev;  noise(k).ac = ac;  noise(k).z = zd;
    noise(k).kurt = kurtosis(zd);

    fprintf('%-5s: %.1f s (%d samples), mean %.2f, std %.3f counts (%.4f deg), p-p %d counts, kurtosis %.2f\n', ...
        R(k).label, noise(k).dur, numel(z), R(k).c0, noise(k).std, noise(k).std*DEG_PER_COUNT, noise(k).p2p, noise(k).kurt);
    fprintf('       drift %.3f counts/s, lag-1 autocorrelation %.3f, Allan dev at 0.25 s = %.3f counts\n', ...
        lin(1), ac(2), interp1(tau, adev, 0.25));
end
fprintf('Effective resolution at 2 kHz: %.1f bits (14 - log2(sigma/q_ideal))\n', ENC_BITS - log2(mean([noise.std])*sqrt(12)));
fprintf('Rest position shift (start -> end of run, at zero airspeed; pivot friction):\n');
for k = 1:numel(R)
    fprintf('  %-5s: start %.1f, end %.1f  (shift %+.1f counts = %+.2f deg)\n', ...
        R(k).label, R(k).c0, R(k).cEnd, R(k).cEnd - R(k).c0, (R(k).cEnd - R(k).c0)*DEG_PER_COUNT);
end


%% Align sensor to pitot
% Trailing 0.25 s mean of counts at each pitot timestamp. A lag scan found
% no measurable delay between the two (< 0.05 s), so none is applied.
T = table();
for k = 1:numel(R)
    cm = interp1(R(k).ts, movmean(R(k).c, [round(WIN_S*R(k).fs)-1 0]), R(k).tg);
    dth = -interp1(R(k).ts, gradient(movmean(R(k).c, round(0.5*R(k).fs)), R(k).ts), R(k).tg);
    ok = ~isnan(cm);
    T = [T; table(R(k).tg(ok) - R(k).ts(1), R(k).q(ok), R(k).vgt(ok), R(k).c0 - cm(ok), dth(ok), ...
        repmat(k, nnz(ok), 1), 'VariableNames', {'t', 'q', 'v', 'th', 'dth', 'run'})]; %#ok<AGROW>
end
fitMask = T.q >= Q_MIN_FIT;
vOf = @(q) sqrt(2*max(q, 0)/rho);
radPerCount = 2*pi / ENC_RES;

%% Candidate models
% Each predicts q from theta (counts, relative to that run's zero tare);
% v = sqrt(2*max(q,0)/rho).
models = {
    'sqrt: v = k*sqrt(theta)',                 @(th) th,                         @(b, th) b(1)*th
    'sqrt+offset: v = k*sqrt(theta-theta0)',   @(th) [th, ones(size(th))],       @(b, th) b(1)*th + b(2)
    'quadratic: q = a1*theta + a2*theta^2',    @(th) [th, th.^2],                @(b, th) b(1)*th + b(2)*th.^2
    'quadratic+offset: q = a0+a1*th+a2*th^2',  @(th) [ones(size(th)), th, th.^2], @(b, th) b(1) + b(2)*th + b(3)*th.^2
    'projected area: q = a*theta/cos^2(theta)', @(th) th./cos(th*radPerCount).^2, @(b, th) b(1)*th./cos(th*radPerCount).^2
    };
nM = size(models, 1);
iSqrt = 1;  iQuad = 3;  iArea = 5;

fprintf('\n===== MODEL COMPARISON (fit on q >= %g Pa, i.e. v >= %.1f m/s; %d pitot samples) =====\n', ...
    Q_MIN_FIT, vOf(Q_MIN_FIT), nnz(fitMask));
fprintf('%-44s %4s %8s %8s %10s %9s %8s\n', 'Model', 'npar', 'RMSE_v', 'RMSE_q', 'CV RMSE_v', 'bias_v', 'BIC');
fprintf('%-44s %4s %8s %8s %10s %9s\n', '', '', '(m/s)', '(Pa)', '(m/s)', '(m/s)');
res = struct([]);
nFit = nnz(fitMask);
for i = 1:nM
    X = models{i, 2};  f = models{i, 3};
    b = X(T.th(fitMask)) \ T.q(fitMask);
    qhat = f(b, T.th);
    % Leave-one-run-out: train on one run, test on the other
    cvErr = [];
    for k = 1:numel(R)
        tr = fitMask & T.run ~= k;  te = fitMask & T.run == k;
        bk = X(T.th(tr)) \ T.q(tr);
        cvErr = [cvErr; vOf(f(bk, T.th(te))) - T.v(te)]; %#ok<AGROW>
    end
    ev = vOf(qhat(fitMask)) - T.v(fitMask);
    eq = qhat(fitMask) - T.q(fitMask);
    res(i).name = models{i, 1};  res(i).b = b;  res(i).f = f;  res(i).qhat = qhat;
    res(i).rmse_v = rms(ev);  res(i).rmse_q = rms(eq);
    res(i).cv = rms(cvErr);   res(i).bias = mean(ev);
    res(i).bic = nFit*log(mean(eq.^2)) + numel(b)*log(nFit);
    fprintf('%-44s %4d %8.3f %8.3f %10.3f %+9.3f %8.1f\n', res(i).name, numel(b), res(i).rmse_v, res(i).rmse_q, ...
        res(i).cv, res(i).bias, res(i).bic);
end
fprintf('(BIC: lower is better; a difference > 10 is conventionally decisive.)\n');

%% The two reported models
a  = res(iSqrt).b;                         % sqrt model
a1 = res(iQuad).b(1);  a2 = res(iQuad).b(2);
covSqrt = ols_cov(T.th(fitMask), T.q(fitMask), a);
covQuad = ols_cov([T.th(fitMask), T.th(fitMask).^2], T.q(fitMask), res(iQuad).b);

M = struct([]);
M(1).name = 'sqrt';       M(1).idx = iSqrt;
M(1).qOf  = @(th) a*th;
M(1).dq   = @(th) a*ones(size(th));
M(1).thOf = @(v) 0.5*rho*v.^2 / a;
M(2).name = 'quadratic';  M(2).idx = iQuad;
M(2).qOf  = @(th) a1*th + a2*th.^2;
M(2).dq   = @(th) a1 + 2*a2*th;
M(2).thOf = @(v) (-a1 + sqrt(a1^2 + 2*a2*rho*v.^2)) / (2*a2);
for j = 1:2, M(j).vOf = @(th) vOf(M(j).qOf(th)); end

fprintf('\n===== REPORTED MODELS  (theta = c0 - count, c0 = zero-airspeed tare) =====\n');
fprintf('[sqrt]       q = a*theta                v = k*sqrt(theta)\n');
fprintf('             a = %.5f +/- %.5f Pa/count\n', a, sqrt(covSqrt));
fprintf('             k = sqrt(2a/rho) = %.4f m/s per sqrt(count)  (%.3f m/s per sqrt(deg))\n', ...
    sqrt(2*a/rho), sqrt(2*a/rho)/sqrt(DEG_PER_COUNT));
fprintf('[quadratic]  q = a1*theta + a2*theta^2  v = k*sqrt(theta)*sqrt(1 + beta*theta)\n');
fprintf('             a1 = %.5f +/- %.5f Pa/count,  a2 = %.3e +/- %.1e Pa/count^2\n', a1, sqrt(covQuad(1,1)), a2, sqrt(covQuad(2,2)));
fprintf('             k = sqrt(2*a1/rho) = %.4f m/s per sqrt(count),  beta = a2/a1 = %.3e /count\n', sqrt(2*a1/rho), a2/a1);
fprintf('             sqrt-law correction factor at theta = 250 / 500 / 950 counts: x%.3f / x%.3f / x%.3f\n', sqrt(1 + a2/a1*[250 500 950]));

fprintf('\n%-12s %8s %8s %10s   %s\n', 'Model', 'RMSE_v', 'RMSE_q', 'CV RMSE_v', 'RMSE_v by speed band (m/s): 2-4 | 4-6 | 6-8 | 8-10.5');
vb = [2 4 6 8 10.5];
for j = 1:2
    r = res(M(j).idx);  ev = vOf(r.qhat) - T.v;
    fprintf('%-12s %8.3f %8.3f %10.3f   ', M(j).name, r.rmse_v, r.rmse_q, r.cv);
    for bb = 1:numel(vb)-1
        mb = fitMask & T.v >= vb(bb) & T.v < vb(bb+1);
        fprintf(' %5.3f |', rms(ev(mb)));
    end
    fprintf('\n');
    for k = 1:numel(R)
        fprintf('%12s %-5s run: mean q residual %+.2f Pa, v RMSE %.3f m/s\n', '', R(k).label, ...
            -mean(r.qhat(fitMask & T.run == k) - T.q(fitMask & T.run == k)), rms(ev(fitMask & T.run == k)));
    end
end

%% Justification of the quadratic term
% The question is whether the curvature is a real property of the sensor or
% an artefact of the known pivot friction (zero shift between runs and
% rising/falling hysteresis) or of one particular run.
fprintf('\n===== JUSTIFICATION OF THE QUADRATIC TERM =====\n');
thF = T.th(fitMask);  qF = T.q(fitMask);  runF = T.run(fitMask);
dirn = latch_direction(T.dth(fitMask), runF, 5);

fprintf('1) Fitted to each run independently:\n');
for k = 1:numel(R)
    mk = runF == k;  X = [thF(mk), thF(mk).^2];
    bk = X \ qF(mk);  ck = ols_cov(X, qF(mk), bk);
    fprintf('     %-5s: a2 = %.3e +/- %.1e Pa/count^2  (t = %.0f)\n', R(k).label, bk(2), sqrt(ck(2,2)), bk(2)/sqrt(ck(2,2)));
end

fprintf('2) With friction terms as nuisance parameters (per-run zero offset + rising/falling hysteresis):\n');
Xn = [double(runF == 1), double(runF == 2), dirn];
bq = [thF, thF.^2, Xn] \ qF;  cq = ols_cov([thF, thF.^2, Xn], qF, bq);
bl = [thF, Xn] \ qF;
rq = qF - [thF, thF.^2, Xn]*bq;  rl = qF - [thF, Xn]*bl;
fprintf('     quadratic: a2 = %.3e +/- %.1e (t = %.0f), RMSE_q %.2f Pa, hysteresis +/- %.2f Pa\n', ...
    bq(2), sqrt(cq(2,2)), bq(2)/sqrt(cq(2,2)), rms(rq), bq(5));
fprintf('     sqrt     : RMSE_q %.2f Pa, needs offsets of %.1f / %.1f Pa (= deadband of %.0f / %.0f counts)\n', ...
    rms(rl), bl(2), bl(3), -bl(2)/bl(1), -bl(3)/bl(1));

fprintf('3) Power law q = c*theta^n fitted per run: ');
for k = 1:numel(R)
    pk = polyfit(log(thF(runF == k)), log(qF(runF == k)), 1);
    fprintf('%s n = %.2f   ', R(k).label, pk(1));
end
fprintf('(sqrt law requires n = 1)\n');

fprintf('4) Sensitivity q/theta (Pa/count) by deflection, constant for the sqrt law:\n');
thBins = 100:100:1000;
sens = nan(numel(thBins)-1, numel(R));
for k = 1:numel(R)
    fprintf('     %-5s:', R(k).label);
    for bb = 1:numel(thBins)-1
        mb = runF == k & thF >= thBins(bb) & thF < thBins(bb+1);
        if nnz(mb) >= 5, sens(bb, k) = median(qF(mb) ./ thF(mb)); end
        fprintf(' %6.4f', sens(bb, k));
    end
    fprintf('\n');
end
fprintf('     bins (counts): %s\n', sprintf(' %4d-%-4d', [thBins(1:end-1); thBins(2:end)]));

fprintf('5) Physical candidates for the curvature:\n');
fprintf('     projected area q = a*theta/cos^2(theta): RMSE_q %.2f Pa, gives only x%.2f sensitivity rise at 950 counts (data: x%.2f)\n', ...
    res(iArea).rmse_q, 1/cos(950*radPerCount)^2, 1 + a2/a1*950);
phis = 0:0.5:40;  rr = zeros(size(phis));
for ii = 1:numel(phis)
    gph = thF ./ cos(thF*radPerCount + deg2rad(phis(ii))).^2;  rr(ii) = rms(qF - gph*(gph \ qF));
end
[rrMin, ii] = min(rr);
fprintf('     same with a rest-angle offset phi (vane not normal to flow at rest): best phi = %.1f deg, RMSE_q %.2f Pa\n', phis(ii), rrMin);

% Binned residual curves for the justification figure
resBins = 0:50:1000;  resCurve = nan(numel(resBins)-1, 2, numel(R));
for j = 1:2
    rj = T.q - res(M(j).idx).qhat;
    for k = 1:numel(R)
        for bb = 1:numel(resBins)-1
            mb = fitMask & T.run == k & T.th >= resBins(bb) & T.th < resBins(bb+1);
            if nnz(mb) >= 4, resCurve(bb, j, k) = mean(rj(mb)); end
        end
    end
end

%% Noise through the board pipeline
% 2 kHz -> 10-sample block mean -> 200 Hz -> N-sample moving average.
% sigma_v = dv/dtheta * sigma_theta, with dv/dtheta = q'(theta) / (rho*v).
fprintf('\n===== NOISE THROUGH BOARD PIPELINE (zero-flow encoder noise, ramp run) =====\n');
zAll = noise(1).z;
zDec = block_mean(zAll, DECIM);
stages = {'2 kHz raw', zAll, 0};
for N = OUT_AVG
    zm = movmean(zDec, [N-1 0]);
    stages(end+1, :) = {sprintf('200 Hz, %2d-tap (%2.0f ms)', N, 5*N), zm(N:end), ...
        1e3*((DECIM-1)/2/R(1).fs + (N-1)/2/(R(1).fs/DECIM))}; %#ok<SAGROW>
end
sigTh = cellfun(@std, stages(:, 2));
vq = [2 3 5 8 10];
for j = 1:2
    fprintf('[%s] sigma_v (m/s):\n', M(j).name);
    fprintf('  %-26s %10s %9s', 'Stage', 'sig_theta', 'delay ms'); fprintf('  v=%4.1f', vq); fprintf('\n');
    dvdth = M(j).dq(M(j).thOf(vq)) ./ (rho*vq);
    for st = 1:size(stages, 1)
        fprintf('  %-26s %10.3f %9.1f', stages{st, 1}, sigTh(st), stages{st, 3}); fprintf('  %6.4f', dvdth*sigTh(st)); fprintf('\n');
    end
    fprintf('  One-count step at v = 2 / 5 / 10 m/s = %.3f / %.3f / %.3f m/s\n', M(j).dq(M(j).thOf([2 5 10])) ./ (rho*[2 5 10]));
end
fprintf('(white noise -> sigma scales as 1/sqrt(samples averaged): %.3f/sqrt(100) = %.3f expected at 50 ms)\n', ...
    sigTh(1), sigTh(1)/10);
fprintf('Delay = group delay of block mean + moving average.\n');

%% Spectrum under flow
% Checks for tonal content (vortex shedding, fan blade-pass, vane resonance)
% that a notch or designed FIR could target, and for aliasing through the
% on-board block mean. Zero-flow segment of the step run is the reference.
fprintf('\n===== SPECTRUM UNDER FLOW (2 kHz counts, linear detrend per segment) =====\n');
segs  = [2, 0, (R(2).iz-1)/R(2).fs; PLATEAUS];
bands = [0.25 2; 2 10; 10 30; 30 inf];
nfftLo = 8192;  nfftHi = 2048;           % 0.25 Hz and 1 Hz resolution
fcMA = 0.443 * R(1).fs / (DECIM*OUT_AVG(end));
dvQuad = @(v) M(2).dq(M(2).thOf(v)) ./ (rho*v);
fprintf('%-5s %-11s %6s %7s | %-31s | %-13s | %s\n', 'Run', 'Window (s)', 'v', 'pitot', ...
    'Variance by band (counts^2)', 'Largest peak', '50 ms output wobble');
fprintf('%-5s %-11s %6s %7s | %7s %7s %7s %7s | %-13s | %s\n', '', '', '(m/s)', 'std Pa', ...
    '0.25-2', '2-10', '10-30', '>30Hz', '>30 Hz', '(counts, m/s)');
spec = struct([]);
for j = 1:size(segs, 1)
    k = segs(j,1);  t = R(k).ts - R(k).ts(1);
    m = t >= segs(j,2) & t < segs(j,3);
    x = detrend(R(k).c(m), 1);
    [pLo, fLo] = pwelch(x, hann(min(nfftLo, floor(numel(x)/2))), [], nfftLo, R(k).fs);
    [pHi, fHi] = pwelch(x, hann(nfftHi), nfftHi/2, nfftHi, R(k).fs);
    bv = arrayfun(@(i) sum(pLo(fLo >= bands(i,1) & fLo < bands(i,2))) * (fLo(2) - fLo(1)), 1:size(bands, 1));

    % Tone check: largest bin above the local median floor
    ex = 10*log10(pHi ./ movmedian(pHi, 51));  hi = fHi > 30;  fh = fHi(hi);
    [exMax, iMax] = max(ex(hi));

    % Fluctuation left in the board output (10-sample block mean + 10-tap MA)
    cOut = movmean(block_mean(R(k).c(m), DECIM), [OUT_AVG(end)-1 0]);
    sOut = std(detrend(cOut(OUT_AVG(end):end), 1));
    tg = R(k).tg - R(k).ts(1);  mg = tg >= segs(j,2) & tg < segs(j,3);
    v = vOf(median(R(k).q(mg)));

    spec(j).label = sprintf('%s %.1f m/s', R(k).label, v);
    spec(j).fLo = fLo;  spec(j).pLo = pLo;  spec(j).fHi = fHi;  spec(j).pHi = pHi;  spec(j).v = v;
    if v > 0, svTxt = sprintf('%.2f / %.4f', sOut, sOut*dvQuad(v)); else, svTxt = sprintf('%.2f / -', sOut); end
    fprintf('%-5s %5.1f-%-5.1f %6.1f %7.2f | %7.2f %7.2f %7.2f %7.2f | +%4.1f dB %4.0f Hz | %s\n', R(k).label, ...
        segs(j,2), segs(j,3), v, std(R(k).q(mg)), bv, exMax, fh(iMax), svTxt);
end
fprintf('Largest peak: highest 1 Hz bin above the local median floor. Compare with the zero-flow row:\n');
fprintf('  a real tone would stand well clear of it and repeat across segments.\n');
fprintf('50 ms output wobble: std of the 200 Hz, 10-tap output (detrended); m/s via the quadratic model.\n');
fprintf('  Encoder floor alone is %.2f counts. Pipeline -3 dB point = %.1f Hz.\n', sigTh(end), fcMA);

%% Save parameters
params = struct( ...
    'sqrt',      struct('a_Pa_per_count', a, 'k_mps_per_sqrtcount', sqrt(2*a/rho)), ...
    'quadratic', struct('a1_Pa_per_count', a1, 'a2_Pa_per_count2', a2, ...
                        'k_mps_per_sqrtcount', sqrt(2*a1/rho), 'beta_per_count', a2/a1), ...
    'rho', rho, 'enc_bits', ENC_BITS, 'noise_std_counts', mean([noise.std]), ...
    'note', 'theta = c0 - count; q = a*theta (sqrt) or a1*theta + a2*theta^2 (quadratic); v = sqrt(2*max(q,0)/rho)');
save(fullfile(rootDir, 'airspeed_model_params.mat'), 'params');

%% Figures (saved to file only)
mcol = [0.85 0.33 0.10; 0.49 0.18 0.56];     % sqrt, quadratic
rcol = lines(numel(R));

% 1. Noise characterisation
fig = figure('Position', [50 50 1200 800]);
tiledlayout(2, 2, 'Padding', 'compact');
nexttile; hold on;
for k = 1:numel(noise)
    zz = noise(k).z;
    histogram(zz, (floor(min(zz))-0.5):(ceil(max(zz))+0.5), 'Normalization', 'pdf', 'DisplayName', R(k).label);
end
xg = linspace(-8, 8, 200); plot(xg, normpdf(xg, 0, noise(1).std), 'k', 'LineWidth', 1.5, 'DisplayName', 'Gaussian fit (ramp)');
xlabel('Count deviation from mean'); ylabel('PDF'); title('Zero-flow count distribution'); legend; grid on;
nexttile; hold on;
for k = 1:numel(noise), plot(noise(k).f, 10*log10(noise(k).pxx), 'DisplayName', R(k).label); end
yline(10*log10((1/12)/(R(1).fs/2)), 'k--', 'Ideal 14-bit quantisation', 'DisplayName', 'Quantisation floor');
xlabel('Frequency (Hz)'); ylabel('PSD (dB counts^2/Hz)'); title('Zero-flow noise PSD (Welch)'); legend; grid on;
nexttile;
for k = 1:numel(noise), loglog(noise(k).tau, noise(k).adev, 'o-', 'DisplayName', R(k).label); hold on; end
loglog(noise(1).tau, noise(1).std ./ sqrt(noise(1).tau*R(1).fs), 'k--', 'DisplayName', 'White noise (\tau^{-1/2})');
xlabel('Averaging time \tau (s)'); ylabel('Allan deviation (counts)'); title('Allan deviation'); legend; grid on;
nexttile; hold on;
for k = 1:numel(noise), stem(0:20, noise(k).ac, 'filled', 'DisplayName', R(k).label); end
xlabel('Lag (samples)'); ylabel('Autocorrelation'); title('Zero-flow autocorrelation'); legend; grid on;
save_fig(fig, figDir, 'noise_characterisation');

% 2. Calibration scatter and fits
fig = figure('Position', [80 80 1200 500]);
tiledlayout(1, 2, 'Padding', 'compact');
thg = linspace(0, max(T.th)*1.05, 300)';
for dom = 1:2
    nexttile; hold on;
    for k = 1:numel(R)
        m = T.run == k;
        if dom == 1, y = T.q(m); else, y = T.v(m); end
        scatter(T.th(m), y, 10, rcol(k,:), 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', R(k).label);
    end
    for j = 1:2
        y = M(j).qOf(thg);  if dom == 2, y = vOf(y); end
        plot(thg, y, 'Color', mcol(j,:), 'LineWidth', 1.8, 'DisplayName', res(M(j).idx).name);
    end
    xlabel('\theta = c_0 - count (counts)'); grid on;
    if dom == 1
        yline(Q_MIN_FIT, ':', 'fit threshold', 'HandleVisibility', 'off');
        ylabel('Pitot dynamic pressure q (Pa)'); title('Pressure domain'); legend('Location', 'northwest');
    else
        ylabel('Pitot airspeed (m/s)'); title('Velocity domain'); legend('Location', 'southeast');
    end
end
save_fig(fig, figDir, 'calibration_fit');

% 3. Justification of the quadratic term
fig = figure('Position', [100 100 1200 450]);
tiledlayout(1, 2, 'Padding', 'compact');
nexttile; hold on;
thc = (thBins(1:end-1) + thBins(2:end)) / 2;
for k = 1:numel(R)
    mk = runF == k;
    scatter(thF(mk), qF(mk) ./ thF(mk), 8, rcol(k,:), 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
    plot(thc, sens(:, k), 'o-', 'Color', rcol(k,:), 'MarkerFaceColor', rcol(k,:), 'LineWidth', 1.3, 'DisplayName', [R(k).label ' (binned median)']);
end
plot(thg, a*ones(size(thg)), 'Color', mcol(1,:), 'LineWidth', 1.8, 'DisplayName', 'sqrt model (constant)');
plot(thg, a1 + a2*thg, 'Color', mcol(2,:), 'LineWidth', 1.8, 'DisplayName', 'quadratic model (a_1 + a_2\theta)');
xlim([0 1000]); ylim([0 0.12]);
xlabel('\theta (counts)'); ylabel('q / \theta (Pa/count)'); title('Sensitivity rises with deflection'); legend('Location', 'northwest'); grid on;
nexttile; hold on;
rbc = (resBins(1:end-1) + resBins(2:end)) / 2;
ls = {'-', '--'};
for j = 1:2
    for k = 1:numel(R)
        plot(rbc, resCurve(:, j, k), ls{k}, 'Color', mcol(j,:), 'Marker', 'o', 'MarkerSize', 4, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s, %s run', M(j).name, R(k).label));
    end
end
yline(0, 'k', 'HandleVisibility', 'off'); xlabel('\theta (counts)'); ylabel('Mean q residual, pitot - model (Pa)');
title('Binned residuals: sqrt model bows the same way in both runs'); legend('Location', 'northwest'); grid on;
save_fig(fig, figDir, 'model_justification');

% 4. Residuals in velocity
fig = figure('Position', [110 110 1200 450]);
tiledlayout(1, 2, 'Padding', 'compact');
vv = linspace(vOf(Q_MIN_FIT), max(T.v), 100);  sq = 1/sqrt(12);
for j = 1:2
    nexttile; hold on;
    ev = vOf(res(M(j).idx).qhat) - T.v;
    for k = 1:numel(R)
        m = fitMask & T.run == k;
        scatter(T.v(m), ev(m), 10, rcol(k,:), 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', R(k).label);
    end
    plot(vv, sq./(rho*vv), 'k--', 'DisplayName', '\pm1\sigma pitot quantisation'); plot(vv, -sq./(rho*vv), 'k--', 'HandleVisibility', 'off');
    yline(0, 'k', 'HandleVisibility', 'off'); ylim([-1 1.2]);
    xlabel('Pitot airspeed (m/s)'); ylabel('Model - pitot (m/s)');
    title(sprintf('%s: RMSE %.3f m/s (CV %.3f)', res(M(j).idx).name, res(M(j).idx).rmse_v, res(M(j).idx).cv)); legend; grid on;
end
save_fig(fig, figDir, 'residuals');

% 5. Time series through the board pipeline (block mean + 10-tap at 200 Hz)
fig = figure('Position', [140 140 1200 700]);
tiledlayout(numel(R), 1, 'Padding', 'compact');
for k = 1:numel(R)
    nexttile; hold on;
    cOut = movmean(block_mean(R(k).c, DECIM), [OUT_AVG(end)-1 0]);
    tOut = block_mean(R(k).ts, DECIM) - R(k).ts(1);
    thO  = R(k).c0 - cOut;
    for j = 1:2, plot(tOut, M(j).vOf(thO), 'Color', mcol(j,:), 'LineWidth', 1.1, 'DisplayName', res(M(j).idx).name); end
    stairs(R(k).tg - R(k).ts(1), R(k).vgt, 'k', 'DisplayName', 'Pitot');
    ylabel('Airspeed (m/s)'); title(sprintf('%s run (200 Hz, 50 ms moving average)', R(k).label)); legend('Location', 'northwest'); grid on;
end
xlabel('Time (s)');
save_fig(fig, figDir, 'timeseries_comparison');

% 6. Propagated airspeed noise through the pipeline
fig = figure('Position', [170 170 1100 450]);
tiledlayout(1, 2, 'Padding', 'compact');
vv = linspace(1, 11, 200);
for j = 1:2
    nexttile; hold on;
    dvdth = M(j).dq(M(j).thOf(vv)) ./ (rho*vv);
    for st = 1:size(stages, 1)
        plot(vv, dvdth*sigTh(st), 'LineWidth', 1.3, 'DisplayName', stages{st, 1});
    end
    set(gca, 'YScale', 'log'); ylim([1e-4 0.1]);
    xlabel('Airspeed (m/s)'); ylabel('\sigma_v from encoder noise (m/s)');
    title(sprintf('%s model', M(j).name)); legend; grid on;
end
save_fig(fig, figDir, 'airspeed_noise');

% 7. Spectra on steady segments
fig = figure('Position', [200 200 1200 800]);
tiledlayout(2, 1, 'Padding', 'compact');
scol = turbo(numel(spec) + 1);
nexttile; hold on;
for j = 1:numel(spec)
    plot(spec(j).fHi, 10*log10(movmean(spec(j).pHi, 3)), 'Color', scol(j,:), 'DisplayName', spec(j).label);
end
yline(10*log10((1/12)/(R(1).fs/2)), 'k--', 'Ideal 14-bit quantisation', 'HandleVisibility', 'off');
xlim([0 R(1).fs/2]); ylim([-45 5]);
xlabel('Frequency (Hz)'); ylabel('PSD (dB counts^2/Hz)'); title('Full band: flat above ~30 Hz at every speed, no tones');
legend('NumColumns', 4, 'Location', 'north'); grid on;
nexttile; hold on;
for j = 1:numel(spec)
    plot(spec(j).fLo, 10*log10(spec(j).pLo), 'Color', scol(j,:), 'LineWidth', 1.1, 'DisplayName', spec(j).label);
end
xline(fcMA, 'k--', sprintf('50 ms output -3 dB (%.1f Hz)', fcMA), 'HandleVisibility', 'off');
set(gca, 'XScale', 'log'); xlim([0.25 200]); ylim([-35 15]);
xlabel('Frequency (Hz)'); ylabel('PSD (dB counts^2/Hz)');
title('Low band (0.25 Hz resolution): broadband flow content below ~20 Hz, growing with speed'); grid on;
save_fig(fig, figDir, 'flow_spectra');

% 8. Spectrograms of both runs with pitot speed
fig = figure('Position', [230 230 1200 800]);
tiledlayout(numel(R), 1, 'Padding', 'compact');
for k = 1:numel(R)
    nexttile;
    x = R(k).c - movmean(R(k).c, round(0.5*R(k).fs));
    [S, F, Ts] = spectrogram(x, hann(2048), 1536, 2048, R(k).fs);
    yyaxis left;
    imagesc(Ts, F, 10*log10(abs(S).^2 / (R(k).fs*sum(hann(2048).^2)))); axis xy;
    clim([-40 -10]); colormap(parula); ylabel('Frequency (Hz)');
    cb = colorbar; cb.Label.String = 'dB counts^2/Hz';
    yyaxis right;
    stairs(R(k).tg - R(k).ts(1), R(k).vgt, 'w', 'LineWidth', 1.3);
    ylim([0 20]); ylabel('Pitot airspeed (m/s)');
    title(sprintf('%s run spectrogram (counts, 0.5 s high-passed)', R(k).label));
end
xlabel('Time (s)');
save_fig(fig, figDir, 'spectrogram');

fprintf('\nFigures written to %s\n', figDir);

%% Local functions
function [tau, adev] = allan_dev(x, fs)
% Non-overlapping Allan deviation over octave-spaced cluster sizes.
    N = numel(x);
    m = unique(round(2.^(0:0.5:floor(log2(N/8)))));
    adev = zeros(size(m));
    for i = 1:numel(m)
        nb = floor(N / m(i));
        y  = mean(reshape(x(1:nb*m(i)), m(i), nb), 1);
        adev(i) = sqrt(0.5 * mean(diff(y).^2));
    end
    tau = m(:) / fs;  adev = adev(:);
end

function y = block_mean(x, n)
% Non-overlapping n-sample means, as done on the board before serial output.
    nb = floor(numel(x) / n);
    y  = mean(reshape(x(1:nb*n), n, nb), 1)';
end

function C = ols_cov(X, y, b)
% Parameter covariance for an ordinary least-squares fit.
    r = y - X*b;
    C = (r'*r) / (numel(y) - numel(b)) * inv(X'*X); %#ok<MINV>
end

function d = latch_direction(dth, run, thresh)
% +1 while theta has most recently been rising, -1 while falling.
    d = ones(size(dth));  last = 1;
    for i = 1:numel(dth)
        if i > 1 && run(i) ~= run(i-1), last = 1; end
        if dth(i) > thresh, last = 1; elseif dth(i) < -thresh, last = -1; end
        d(i) = last;
    end
end

function save_fig(fig, figDir, name)
    exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 150);
    close(fig);
end
