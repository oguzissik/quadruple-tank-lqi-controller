```bash
cat > README.md << 'EOF'
# Quadruple Tank — LQI Control System

Model-based optimal control design and real hardware deployment for a nonlinear coupled quadruple tank system using LQI (Linear Quadratic Integral) control.

## System

Quanser Coupled Tanks — 4 tanks, 2 pump inputs, 4 water level outputs. The tanks are physically coupled: changing one pump voltage affects multiple tank levels simultaneously, making this a challenging MIMO control problem.

## Method

**Grey-box identification** — pump flow and valve opening parameters identified from real experimental data via constrained optimization.

**Linearization** — nonlinear dynamics linearized around an operating equilibrium point using Taylor expansion, producing state-space matrices A, B, C, D.

**LQI design** — augmented system with integral states constructed, optimal gains Kx and Ke computed via LQR on the augmented system.

**Discrete-time implementation** — controller discretized at τs = 0.01s for real hardware deployment.

## Files

- `main_validation.m` — grey-box model validation against experimental data
- `main_control.m` — equilibrium computation, linearization, LQI design, discretization
- `LQI_Controller.slx` — full closed-loop Simulink model with nonlinear tank plant
- `LQI_params.mat` — all necessary parameters for hardware deployment

## Requirements

- MATLAB R2021b or later
- Simulink
- Curve Fitting Toolbox (required to load pump parameter files)
- Optimization Toolbox (required for equilibrium computation via fmincon)

## How to Run

1. Run `main_validation.m` to validate the grey-box model
2. Run `main_control.m` to compute all control parameters
3. Open and simulate `LQI_Controller.slx`

## Results

Controller achieves zero steady-state error with smooth step tracking on both simulation and real Quanser hardware.
EOF

git add README.md
git commit -m "Add README"
git push
```

Hepsini kopyala, terminale yapıştır — tek seferde çalışır.
