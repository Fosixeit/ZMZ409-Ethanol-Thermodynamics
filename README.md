# ZMZ-409 Ethanol Thermodynamics Model

## Overview
This repository contains the open-source Python code for the thermodynamic 0D/1D mathematical model developed to simulate the working cycle of the **ZMZ-409 spark-ignition engine** running on ethanol-gasoline blends (E10, E20, E30).

This code is provided as Supplementary Material to ensure full reproducibility of the research results presented in our paper.

## Features
* **Kinematic Module.**
* Calculates instantaneous cylinder volume based on engine geometry.
* **Combustion Module (Wiebe Function).**
* Adjusts combustion duration and heat release rates for different ethanol mass fractions.
* **Heat Transfer.**
* Implements the Woschni heat transfer model.
* **Interactive Visualization.**
* Features a real-time UI with sliders to test the impact of *Ethanol Fraction* and *Spark Advance (Ignition Timing)* on the $P-\theta$ and $P-V$ diagrams.

## Dependencies
To run the simulation, you need Python 3.x and the following scientific libraries:

```bash
pip install numpy scipy matplotlib
```

## How to Run
Simply execute the main Python script. An interactive dashboard will open automatically.

```bash
python Calculating.py
```
