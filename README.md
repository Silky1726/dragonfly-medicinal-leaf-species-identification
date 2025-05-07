# dragonfly-medicinal-leaf-species-identification
Medicinal leaf species identification using metaheuristic algorithms.
# Medicinal Plant Leaf Classification Framework 🌿🧠
**"Levy’s Flight-based Dragonfly Optimized Machine Learning Framework to Identify Medicinal Plants through Leaf Images"**

## 🔍 Overview

This framework identifies medicinal plant species based on leaf images using:
- K-means clustering, optimized by firefly algorithm for Segmentation of ROI from background.
- BRISK for feature extraction
- Modified Dragonfly algorithm with Levy's flight for feature selection
- MLNN (Multilayer Neural Network) for classification

## 🗂 Structure
- `Segmentation.m`: Image preprocessing and segmentation using K-Means + Firefly
- `FeatureExtraction_BRISK.m`: BRISK keypoint detection and descriptor generation
- `FeatureSelection_Dragonfly.m`: Custom Dragonfly algorithm for feature selection
- `Classifier_Training.m`: MLNN model training and evaluation
- `Dataset/`: Contains sample leaf images (PlantVillage + Custom Medicinal Dataset)
- `model.mat`: Trained model file
- `GUI_Frontend.m`: MATLAB GUI interface

## 🖥 System Requirements

- MATLAB R2021b or higher
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

## ⚙️ How to Run

1. Clone this repo:
    ```bash
    git clone https://github.com/YOUR_USERNAME/Medicinal-Leaf-Classification.git
    ```
2. Open MATLAB and add the project folder to the path.
3. Run `GUI_Frontend.m` to launch the GUI or execute individual `.m` files for testing each component.

## 🔁 Reproducibility Settings

- **Segmentation:** K-Means (k=2), Firefly Algorithm (population=50, α=0.5, β=0.2, γ=1.0)
- **Feature Extraction:** BRISK (OpenCV default config)
- **Feature Selection:** Dragonfly (population=30, iterations=100, weights for cohesion/separation/alignment=0.2)
- **Classifier:** MLNN (2 hidden layers with 100 and 50 neurons, ReLU activation)
- **System:** Intel i7, 16GB RAM, MATLAB R2021b

## 📊 Performance

- **PlantVillage Dataset Accuracy:** 88.6%
- **Custom Medicinal Leaf Dataset Accuracy:** 92.4%
- Outperforms or matches CNN models with lower computational requirements

## 📜 Citation

If you use this framework, please cite our paper:

## 📫 Contact

- **Silky Sachar** – silky.sachar@gmail.com  
- **Anuj Kumar** – [Your Email Here]

---

