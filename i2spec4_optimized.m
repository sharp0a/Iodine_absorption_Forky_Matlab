function results = i2spec4_optimized(temp, pres, rlen, wnstart, wnstop, wnstep, ...
                                     line_data, aets_option)
% I2SPEC4_OPTIMIZED - Iodine absorption spectrum model
% Based on: Forkey et al., Appl. Opt. 36(27):6729-6738 (1997)
%
% PHYSICS:
%   Models I₂ B³Π(0u+) ← X¹Σg+ transitions near 532nm with hyperfine structure
%   - Nuclear electric quadrupole (NEQ) splitting [Eq. 3]
%   - Magnetic spin-rotation splitting [Eq. 5]
%   - Result: ~21 hyperfine components (odd J) or ~15 (even J)
%
% INPUTS:
%   temp       - Temperature [K] (typical: 300-400K)
%   pres       - Pressure [Torr] (keep <3 for accuracy)
%   rlen       - Length [cm] (typical: 10-127cm)
%   wnstart    - Starting wavenumber [cm⁻¹]
%   wnstop     - Ending wavenumber [cm⁻¹]
%   wnstep     - Resolution [cm⁻¹]
%   line_data  - From i2lines2: .nlines, .energy, .fcf, .ivl, .ivu, .j, 
%                .branch, .deltaeqq, .deltac
%   aets_option- 1=Tellinghuisen(1982), 2=Princeton 0.99×10⁻³⁶ (recommended)
%
% OUTPUTS:
%   results.wavenumber   - Frequency array [cm⁻¹]
%   results.transmission - Transmission [0-1]
%   results.parameters   - Input parameters + AETS

    %% Store parameters
    params = struct('temp', temp, 'pres', pres, 'rlen', rlen, ...
                    'wnstart', wnstart, 'wnstop', wnstop, 'wnstep', wnstep, ...
                    'aets_option', aets_option);
    
    %% Extract line data
    nlines = line_data.nlines;
    centwn = line_data.energy;
    fcf = line_data.fcf;
    ivl = line_data.ivl;
    ivh = line_data.ivu;
    ij = line_data.j;
    deltaeqq = line_data.deltaeqq;
    deltac = line_data.deltac;
    
    fprintf('Lines: %d, FCF range: [%.2e, %.2e]\n', nlines, min(fcf), max(fcf));
    
    %% Calculate Boltzmann populations [molecules/cm³]
    rnum = calc_nums(temp, pres, ivl, ij, nlines);
    
    %% Calculate hyperfine structure
    [centhypwn, widthhypwn] = calc_hyperfine(centwn, deltaeqq, deltac, ij, temp, nlines);
    
    %% Get AETS
    aets = get_aets((wnstart + wnstop)/2, aets_option);
    params.aets = aets;
    fprintf('AETS: %.2e esu²cm²\n', aets);
    
    %% Peak absorption coefficients [Eq. 8]
    centabs = rlen * (rnum ./ sqrt(temp)) * 4.348e24 * aets .* fcf;
    fprintf('Optical depth range: [%.2e, %.2e]\n', min(centabs), max(centabs));
    
    %% Calculate transmission [Eq. 11]
    [wavenumber, trans] = calc_trans(wnstart, wnstop, wnstep, centabs, ...
                                      centhypwn, widthhypwn, ij, nlines, centwn);
    
    %% Package results
    results = struct('wavenumber', wavenumber, 'transmission', trans, 'parameters', params);
end

%% ========================================================================
function rnum = calc_nums(temp, pres, ivl, ij, nlines)
    % Calculate Boltzmann population densities
    % N(v",J") = Nₜₒₜ(2J"+1)/Q × exp[-hc E(v",J")/kT]
    
    % X-state constants from Gerstenkorn & Luc (1986)
    EL = [0, 213.3023, 425.3742, 636.2103, 845.8034, 1054.1456, 1261.2284, ...
          1467.0428, 1671.5796, 1874.8295, 2076.7829, 2277.4300, 2476.7606, ...
          2674.7643, 2871.4302, 3066.7470, 3260.7032, 3453.2871, 3644.4874, 3834.2932];
    
    BL = [0.03731114, 0.03719670, 0.03708161, 0.03696586, 0.03684942, 0.03673226, ...
          0.03661438, 0.03649575, 0.03637635, 0.03625616, 0.03613514, 0.03601327, ...
          0.03589051, 0.03576683, 0.03564218, 0.03551653, 0.03538980, 0.03526195, ...
          0.03513291, 0.03500261];
    
    DL = 1e-8 * [0.45475, 0.45722, 0.45975, 0.46237, 0.46509, 0.46791, 0.47084, ...
                 0.47388, 0.47704, 0.48031, 0.48371, 0.48723, 0.49090, 0.49471, ...
                 0.49866, 0.50275, 0.50694, 0.51121, 0.51549, 0.51970];
    
    HL = -1e-15 * [0.512752, 0.534032, 0.555643, 0.577957, 0.601224, 0.625607, ...
                   0.651217, 0.678149, 0.706515, 0.736482, 0.768306, 0.802367, ...
                   0.839206, 0.879559, 0.924391, 0.974934, 1.032720, 1.099620, ...
                   1.177872, 1.270126];
    
    % Total number density [molecules/cm³]
    rntot = (pres/760) * 1.013e6 / (1.381e-16 * temp);
    
    % Partition function
    ptfn = 0;
    for iv = 0:19
        for j = 0:200
            gamma = j*(j+1);
            energywn = EL(iv+1) + BL(iv+1)*gamma - DL(iv+1)*gamma^2 + HL(iv+1)*gamma^3;
            rtdegen = 2*j+1;
            spdegen = 15 + 6*mod(j,2);  % 15 even, 21 odd (nuclear spin)
            ptfn = ptfn + rtdegen*spdegen*exp(-1.4388*energywn/temp);
        end
    end
    
    % Vectorized population calculation
    gamma = ij.*(ij+1);
    energywn = EL(ivl+1)' + BL(ivl+1)'.*gamma - DL(ivl+1)'.*gamma.^2 + HL(ivl+1)'.*gamma.^3;
    rtdegen = 2*ij+1;
    spdegen = 15 + 6*mod(ij,2);
    rnum = rntot * (rtdegen ./ ptfn) .* exp(-1.4388*energywn/temp);
    
    fprintf('  Nₜₒₜ=%.2e, Q=%.2e, N range=[%.2e,%.2e]\n', rntot, ptfn, min(rnum), max(rnum));
end

%% ========================================================================
function [centhypwn, widthhypwn] = calc_hyperfine(centwn, deltaeqq, deltac, ij, temp, nlines)
    % Calculate hyperfine-split line centers and Doppler widths
    % Eq. 3 (NEQ): Δν̃ = (ΔeQq/80)[...]
    % Eq. 5 (mag): Δν̃ = ΔC[J(M₁+M₂) + ...]
    % Eq. 10: FWHM = 2.7×10⁻⁸ ν̃ √T
    
    centhypwn = zeros(nlines, 6, 6);
    widthhypwn = zeros(nlines, 6, 6);
    
    for i = 1:nlines
        % Pauli exclusion: M₁≠M₂ for even J
        is_even = mod(ij(i), 2) == 0;
        
        for m1 = 1:6
            m2max = m1 - is_even;
            for m2 = 1:m2max
                rm1 = m1 - 3.5;  % Convert to physical M: -2.5 to +2.5
                rm2 = m2 - 3.5;
                rj = ij(i);
                
                if rj ~= 0
                    % Magnetic term
                    delta_mag = deltac(i) * (rj*(rm1+rm2) + ...
                                0.5*(rm1*(rm1+1) + rm2*(rm2+1)) - 8.75);
                    
                    % NEQ term
                    delta_neq = -(deltaeqq(i)/80) * ...
                                (3*(rm1^2 + rm2^2) + ...
                                 (3/rj)*(rm1*(rm1*(rm1+1)-8.25) + rm2*(rm2*(rm2+1)-8.25)) - 17.5);
                    
                    centhypwn(i,m1,m2) = centwn(i) + delta_mag + delta_neq;
                else
                    centhypwn(i,m1,m2) = centwn(i);
                end
                
                % Doppler width
                widthhypwn(i,m1,m2) = 2.700e-8 * centhypwn(i,m1,m2) * sqrt(temp);
            end
        end
    end
end

%% ========================================================================
function aets = get_aets(wnum, option)
    % Average Electronic Transition Strength |⟨μₑ(R)⟩|²
    % Option 1: Tellinghuisen (1982) tabulated
    % Option 2: Princeton 0.99×10⁻³⁶ (7% correction for 532nm)
    
    strength = [0.56, 0.63, 0.67, 0.72, 0.78, 0.83, 0.92, 0.94, ...
                1.06, 1.05, 1.02, 1.06, 1.01, 1.01, 1.01, 1.23, ...
                1.07, 1.27, 0.95];
    
    if option == 1
        wave = 1e7/wnum;  % Convert to wavelength [nm]
        idx = max(45, min(62, round(wave/10)));
        aets = (strength(idx-44) + (wave/10-idx)*(strength(idx-43)-strength(idx-44))) * 1e-36;
    else
        aets = 0.99e-36;  % Princeton value for 532nm
    end
end

%% ========================================================================
function [wavenumber, trans] = calc_trans(wnstart, wnstop, wnstep, centabs, ...
                                          centhypwn, widthhypwn, ij, nlines, centwn)
    % Calculate transmission via Beer's law [Eq. 11]
    % T(ν̃) = exp[-Σᵢ Gᵢ gᵢ(ν̃)]
    % Optimized with 0.2 cm⁻¹ sliding window
    
    npts = round((wnstop-wnstart)/wnstep);
    trans = ones(npts, 1);
    wavenumber = (wnstart + wnstep:wnstep:wnstart + npts*wnstep)';
    
    ilo = 1; ihi = 2;
    
    for n = 1:npts
        wavenum = wavenumber(n);
        
        % Progress indicator
        if mod(n, round(npts/10)) == 0
            fprintf('  %.0f%%\n', 100*n/npts);
        end
        
        % Slide window: only lines within ±0.2 cm⁻¹
        while ilo <= nlines && (wavenum - centwn(ilo)) >= 0.2
            ilo = ilo + 1;
        end
        while ihi < nlines && (centwn(ihi+1) - wavenum) <= 0.2
            ihi = ihi + 1;
        end
        
        % Accumulate absorption from all hyperfine components
        for i = ilo:min(ihi, nlines)
            is_even = mod(ij(i), 2) == 0;
            for m1 = 1:6
                m2max = m1 - is_even;
                for m2 = 1:m2max
                    if widthhypwn(i,m1,m2) > 0
                        % Normalized detuning
                        delta = (wavenum - centhypwn(i,m1,m2)) / widthhypwn(i,m1,m2);
                        % Gaussian absorption
                        trans(n) = trans(n) * exp(-centabs(i) * exp(-delta^2));
                    end
                end
            end
        end
    end
end
