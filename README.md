# Diagnosability for Cyber-Physical Systems

An abstraction-free and automata-based framework to verify diagnosability over cyber-physical systems using hybrid barrier certificates.

## Overview

This project develops a novel approach to verify the diagnosability of cyber-physical systems through hybrid barrier certificates. Given a desired diagnosability property, we construct a (δ,K)-deterministic finite automaton ((δ,K)-DFA) that captures fault occurrences targeted for diagnosis. The verification problem is then converted to a safety verification problem over a product system.

## Features
- **Abstraction-free verification** of diagnosability properties
- **Hybrid barrier certificate** computation using:
  - Sum-of-Squares (SOS) programming
  - Counter-Example Guided Inductive Synthesis (CEGIS)
- **Diagnoser construction** for diagnosable systems
- **DFA-based framework** for fault detection

## Installation
1. Create Python environment
```bash
conda create -n approx-diag python=3.8 
pip install "cvxpy[MOSEK]"
```

2. Install Julia environment from [Julia](https://julialang.org/downloads/) and add the following packages
```julia
add ReachabilityAnalysis Plots LazySets LaTeXStrings
```

## Usage
1. To reproduce the B-HBC in the Running example, type the following commands
```bash
cd running-example\runing-example-B
python main.py
```
Note that in those output equations obtained in the Terminal, $x1$ and $x2$ are used to denote $x$ and $\hat{x}$ in our manuscript, respectively. Moreover, B0 - B5 correspond to the cases where $q=\bar{q}_0$, $q=q_1$, $q=q_2$, $q=q_3$, $q=q_{\text{trap}}$, and $q=\bar{F}$, respectively.

2. To reproduce the $\mathcal{V}$-HBC in the Running example, type the following commands 
```bash
cd running-example\runing-example-V
python main.py
```
Note that in those output equations obtained in the Terminal, $x1$ and $x2$ are used to denote $x$ and $\hat{x}$ in our manuscript, respectively. Moreover, V0 - V4 correspond to the cases where $q=\bar{q}_0$, $q=q_1$, $q=q_2$, $q=q_{\text{trap}}$, and $q=\bar{F}$, respectively.

3. To reproduce the $\mathcal{B}$-HBC of the Case Study in Section~IV, type the following commands
```bash
cd final-case\B-function-3-pkg
python main.py
```

4. To reproduce the $\mathcal{V}$-HBC of the Case Study in Section~IV, type the following commands
```bash
cd final-case\V-function-1
python main.py
```
Note that in those output equations obtained in the Terminal, $x1$, $x2$, $x12$, and $x22$ are used to denote $x_1$, $x_2$ $\hat{x}_1$, $\hat{x}_2$ in our manuscript, respectively.

5. To reproduce the Figure~4 of the Case Study in Section~IV, type the following commands
```bash
cd final-case\diagnoser
julia
include("diagnoser.jl")
```
