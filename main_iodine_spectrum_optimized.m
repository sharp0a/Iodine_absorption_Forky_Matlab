%% MAIN_IODINE_SPECTRUM - I₂ absorption calculator (optimized)
% Based on Forkey, Lempert, Miles (Applied Optics, 1997)
% User inputs: nm (vaccume) | Calculations: wavenumbers
% Based in fortran Forkey code
%

clear; clc;
%%
opts.Resize = 'on';
opts.WindowStyle = 'normal';

%% Choose line data source
source_choice = questdlg(...
    sprintf(['Line Position Data Source:\n\n' ...
             'Generate New: Calculate from spectroscopic constants\n' ...
             '  • Requires fcfioded file\n' ...
             '  • ~1-2 sec per nm\n\n' ...
             'Load from File: Pre-calculated (faster)\n' ...
             '  • Instant\n' ...
             '  • Use provided .txt files']), ...
    'Select Line Data Source', ...
    'Generate New', 'Load from File', 'Cancel');

if isempty(source_choice) || strcmp(source_choice, 'Cancel')
    disp('Cancelled.'); return;
end

if strcmp(source_choice, 'Generate New')
    %% Generate new line data
    prompts = {'Low wavelength [nm]:', 'High wavelength [nm]:'};
    dlg_title = 'Wavelength Range';
    defaults = {'425', '550'};
    
    answers = inputdlg(prompts, dlg_title, 1, defaults, opts);
    if isempty(answers), disp('Cancelled.'); return; end
    
    lambda_low = str2double(answers{1});
    lambda_high = str2double(answers{2});
    
    if lambda_low >= lambda_high
        error('Low wavelength must be < high wavelength');
    end
    
    wnlo = 1e7 / lambda_high;
    wnhi = 1e7 / lambda_low;
    
    fprintf('Range: %.2f-%.2f nm (%.0f-%.0f cm⁻¹)\n', ...
            lambda_low, lambda_high, wnlo, wnhi);
    
    [fcf_filename, fcf_path] = uigetfile('*', 'Select fcfioded file');
    if isequal(fcf_filename, 0), disp('Cancelled.'); return; end
    fcf_file = fullfile(fcf_path, fcf_filename);
    
    fprintf('Generating line positions...\n'); tic;
    lines = i2lines2_optimized(wnlo, wnhi, fcf_file);
    fprintf('Generated %d transitions in %.2f sec\n', lines.nlines, toc);
    
else
    %% Load pre-generated line data
    [line_filename, line_path] = uigetfile('*.txt', ...
        'Select line data file (e.g., 18786to19791m2.txt)');
    if isequal(line_filename, 0), disp('Cancelled.'); return; end
    line_file = fullfile(line_path, line_filename);
    
    fprintf('Loading line data...\n'); tic;
    fid = fopen(line_file, 'r');
    if fid == -1, error('Cannot open: %s', line_file); end
    
    nlines = fscanf(fid, '%d', 1);
    
    % Pre-allocate for speed
    lines.energy = zeros(nlines, 1);
    lines.ivu = zeros(nlines, 1);
    lines.ivl = zeros(nlines, 1);
    lines.j = zeros(nlines, 1);
    lines.branch = char(zeros(nlines, 1));
    lines.fcf = zeros(nlines, 1);
    lines.deltaeqq = zeros(nlines, 1);
    lines.deltac = zeros(nlines, 1);
    
    for i = 1:nlines
        data = textscan(fid, '%f %d %d %d %s %f %f %f', 1);
        lines.energy(i) = data{1};
        lines.ivu(i) = data{2};
        lines.ivl(i) = data{3};
        lines.j(i) = data{4};
        lines.branch(i) = data{5}{1};
        lines.fcf(i) = data{6};
        lines.deltaeqq(i) = data{7};
        lines.deltac(i) = data{8};
    end
    fclose(fid);
    lines.nlines = nlines;
    
    % Infer wavelength range
    lambda_low = 1e7 / max(lines.energy);
    lambda_high = 1e7 / min(lines.energy);
    
    fprintf('Loaded %d transitions (%.2f-%.2f nm) in %.2f sec\n', ...
            nlines, lambda_low, lambda_high, toc);
end

%% Get spectrum parameters
prompts = {
    'Cell temperature [K]:', ...
    'Cell pressure [Torr]:', ...
    'Cell length [cm]:', ...
    'Start wavelength [nm]:', ...
    'End wavelength [nm]:', ...
    'Resolution [nm]:'
};
dlg_title = 'Spectrum Parameters';
defaults = {'323.15', '2.29', '127.0', ...
            sprintf('%.2f', lambda_low), sprintf('%.2f', lambda_high), '0.0001'};

answers = inputdlg(prompts, dlg_title, 1, defaults, opts);
if isempty(answers), disp('Cancelled.'); return; end

temp = str2double(answers{1});
pres = str2double(answers{2});
rlen = str2double(answers{3});
lambda_start = str2double(answers{4});
lambda_stop = str2double(answers{5});
lambda_step = str2double(answers{6});

if lambda_start >= lambda_stop
    error('Start wavelength must be < end wavelength');
end

% Convert to wavenumbers
wnstart = 1e7 / lambda_stop;
wnstop = 1e7 / lambda_start;
wnstep = (wnstop - wnstart) / ((lambda_stop - lambda_start) / lambda_step);

%% AETS selection (Average Electronic Transition Strength)
aets_choice = questdlg(...
    sprintf(['Average Electronic Transition Strength:\n\n' ...
             'Published: Tellinghuisen (1982) wavelength-dependent\n' ...
             'Princeton: 0.99×10⁻³⁶ esu²cm² (Forkey 1997)\n\n' ...
             'Princeton recommended for 532 nm']), ...
    'AETS Selection', ...
    'Published', 'Princeton', 'Princeton');

if isempty(aets_choice), disp('Cancelled.'); return; end

aets_option = strcmp(aets_choice, 'Published') + 2*strcmp(aets_choice, 'Princeton');

%% Calculate spectrum (optimized version)
fprintf('\n=== Calculating transmission spectrum ===\n');
tic;
spectrum = i2spec4_optimized(temp, pres, rlen, wnstart, wnstop, wnstep, ...
                              lines, aets_option);
calc_time = toc;
fprintf('Calculated %d points in %.2f sec (%.1f pts/sec)\n', ...
        length(spectrum.transmission), calc_time, length(spectrum.transmission)/calc_time);

%% Convert to wavelengths
wavelength_nm = 1e7 ./ spectrum.wavenumber;

%% Results summary
fprintf('\n=== Results ===\n');
fprintf('Lines: %d\n', lines.nlines);
fprintf('Range: %.2f - %.2f nm\n', min(wavelength_nm), max(wavelength_nm));
fprintf('Transmission: %.2e (min), %.4f (max)\n', ...
        min(spectrum.transmission), max(spectrum.transmission));

%% Plot spectrum
fig = figure('Name', 'I₂ Transmission Spectrum', 'NumberTitle', 'off', ...
       'Position', [100 100 1200 800]);

% Optical density
optical_density = -log10(spectrum.transmission);
optical_density(isinf(optical_density)) = max(optical_density(~isinf(optical_density))) + 1;

% Top: Optical Density
subplot(2,1,1);
plot(wavelength_nm, optical_density, 'b-', 'LineWidth', 1.2);
xlabel('Wavelength [nm]', 'FontSize', 11);
ylabel('Optical Density', 'FontSize', 11);
title(sprintf('I₂ Cell: T=%.1f K, P=%.2f Torr, L=%.1f cm', temp, pres, rlen), ...
      'FontSize', 12, 'FontWeight', 'bold');
xlim([min(wavelength_nm) max(wavelength_nm)]);
ylim([0 max(optical_density)*1.05]);
grid on; set(gca, 'FontSize', 10);

% Bottom: Transmission (log)
subplot(2,1,2);
semilogy(wavelength_nm, spectrum.transmission, 'r-', 'LineWidth', 1.2);
xlabel('Wavelength [nm]', 'FontSize', 11);
ylabel('Transmission', 'FontSize', 11);
xlim([min(wavelength_nm) max(wavelength_nm)]);
ylim([max(1e-7, min(spectrum.transmission)*0.1) 1]);
grid on; set(gca, 'FontSize', 10);

% Citation
annotation(fig, 'textbox', [0.15 0.01 0.75 0.03], ...
           'String', 'Forkey et al., Appl. Opt. 36(27):6729-6738 (1997) | Optimized version', ...
           'EdgeColor', 'none', 'FontSize', 8, 'FontAngle', 'italic', ...
           'HorizontalAlignment', 'center');

%% Save option
save_choice = questdlg('Save results?', 'Save', 'Yes', 'No', 'Yes');
if strcmp(save_choice, 'Yes')
    [save_filename, save_path] = uiputfile('*.mat', 'Save As');
    if ~isequal(save_filename, 0)
        results.lines = lines;
        results.spectrum = spectrum;
        results.wavelength_nm = wavelength_nm;
        results.parameters.temp = temp;
        results.parameters.pres = pres;
        results.parameters.rlen = rlen;
        results.parameters.aets_option = aets_option;
        results.metadata.date = datetime('now');
        results.metadata.reference = 'Forkey et al., Appl. Opt. 36(27):6729-6738 (1997)';
        results.metadata.version = 'optimized';
        results.metadata.calc_time_sec = calc_time;
        
        save(fullfile(save_path, save_filename), 'results', '-v7.3');
        fprintf('Saved: %s\n', save_filename);
    end
end

fprintf('\nComplete.\n');

%%