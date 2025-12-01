# Multi-Domain Sensor Fusion for Human Detection
A .NET/WPF-Based Real-Time Multi-Sensor Simulation Framework

This repository contains the implementation of a multi-domain sensor fusion system for human presence detection. The system integrates data from several simulated sensor modalities, including thermal imaging, CO₂ concentration monitoring, radar Doppler signatures, microphone audio data, and IMU motion features.  
The software is developed using **C# (.NET 8.0)** with a **WPF graphical interface**.

## 1. Overview
The system generates real-time simulated sensor data and processes it through a dedicated detection pipeline. Each sensor model contributes unique features, which are combined to estimate the probability of human presence.  
The framework is designed for academic research, algorithm development, and testing.

## 2. Sensors Included

### 2.1 Thermal Sensor (MLX90640 Simulation)
- 32×24 temperature grid  
- Heatmap generation  
- Noise and smoothing models  

### 2.2 CO₂ Sensor
- Baseline CO₂ concentration  
- Variation modeling for human presence detection  

### 2.3 Radar Doppler Sensor
- Motion energy simulation  
- Micro-Doppler pattern generation  
- FFT-based frequency-domain processing  

### 2.4 Microphone Sensor
- RMS amplitude extraction  
- Basic audio event simulation  

### 2.5 IMU Sensor
- Acceleration (X/Y/Z)  
- Orientation and tilt variation  
- Motion variance analysis  

## 3. Detection Pipeline
All sensor features are processed in `DetectionPipeline.cs`. The pipeline includes:

- Threshold-based decision rules  
- Confidence scoring  
- Temporal smoothing  
- Multi-sensor logical fusion  

Example simplified detection rule:

```csharp
bool isHuman =
    (thermal && co2) ||
    (thermal && radar) ||
    (co2 && radar && audio) ||
    (thermal && imu);
