# Iodine_absorption_Forky_Matlab
MATLAB port of the molecular iodine absorption spectrum model by J. N. Forkey

#Quick Start Guide
## Files Provided

1. **i2spec4_optimized.m** - Core calculation engine
2. **main_iodine_spectrum_optimized.m** - User interface script


## Installation

```matlab
% Add to MATLAB path
addpath('/path/to/downloaded/files');
```

## Quick Run

```matlab
% Option 1: Interactive mode
main_iodine_spectrum_optimized

% Option 2: Direct calculation
lines = load_or_generate_lines();  % See below
results = i2spec4_optimized(323.15, 2.29, 127, 18786, 18791, 0.0001, lines, 2);

% Plot
plot(1e7./results.wavenumber, results.transmission);
xlabel('Wavelength [nm]'); ylabel('Transmission');
```

## Loading Line Data

### Option A: From file (fast)
```matlab
fid = fopen('18786to19791m2.txt');
nlines = fscanf(fid, '%d', 1);
lines.energy = zeros(nlines,1);
% ... (see main script lines 50-70)
```

### Option B: Generate new (slow but custom range)
```matlab
lines = i2lines2(wnlo, wnhi, 'fcfioded');
```

## Parameters

| Parameter | Typical Range | Units | Notes |
|-----------|--------------|-------|-------|
| Temperature | 300-400 | K | Higher T → broader lines |
| Pressure | 0.5-3.0 | Torr | Keep <3 for accuracy |
| Length | 10-127 | cm | Longer → stronger absorption |
| Resolution | 0.0001-0.001 | nm | Finer → slower |
| AETS option | 2 | - | Use Princeton for 532nm |

## Example: 532nm Notch Filter

```matlab
% Forkey et al. Fig. 1 conditions
temp = 273.15 + 80;  % 80°C
pres = 0.70;         % Torr
rlen = 25.28;        % cm

% Load pre-calculated lines
lines = load('18786to19791m2.txt');

% Calculate spectrum
tic;
results = i2spec4_optimized(temp, pres, rlen, 18787, 18789, 0.001, lines, 2);
fprintf('Computed in %.1f sec\n', toc);

% Plot
figure;
plot(1e7./results.wavenumber, results.transmission);
xlabel('Wavelength [nm]'); ylabel('Transmission');
title(sprintf('I₂ Filter: %.1f°C, %.2f Torr, %.1f cm', temp-273.15, pres, rlen));
grid on;
```

## Speed Tips

1. **Coarser resolution**: 0.001 nm instead of 0.0001 → 10× faster
2. **Narrow wavelength range**: Only compute what you need
3. **Pre-load lines**: Generate once, save, reuse

## Validation

Compare your results to Forkey et al. (1997) Table 2:
- Line 1: τ₀ ≈ 3.1 (80°C, 0.47 Torr, 25.28 cm)
- Line 17: τ₀ ≈ 13.6

## Troubleshooting

**"Not enough memory"**
- Reduce resolution or wavelength range

**"Transmission = 1 everywhere"**
- Check line data loaded correctly
- Verify wavelength/wavenumber range overlap

**"Lines too narrow"**
- Increase resolution (decrease wnstep)

**"Results don't match paper"**
- Use AETS option 2 (Princeton) for 532nm
- Check temperature in Kelvin, not Celsius

## Citation

If using this code in publications:

> Forkey, J. N., Lempert, W. R., & Miles, R. B. (1997). 
> Corrected and calibrated I₂ absorption model at frequency-doubled Nd:YAG laser wavelengths. 
> Applied Optics, 36(27), 6729-6738.
