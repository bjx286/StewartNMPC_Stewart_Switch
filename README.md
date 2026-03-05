# Stewart 平台 NMPC 仿真项目

本项目包含了 Stewart 平台的非线性模型预测控制（NMPC）仿真流程。

## 🧩 核心代码模块

* **`main`**：主流程与作图
* **`buildConfig`**：参数配置
* **`runSimulation`**：ASM + EKF + 仿真循环
* **`SolVeNmpcQp`**：NMPC QP 求解
* **`evaluateMetrics`**：MIFW 指标评估

## 🚀 运行方式

1. 打开 MATLAB，将当前工作路径切换到以下目录：
   ```text
   d:\TEST\MatLab\AI_generate\Stewart20191126_NMPC_Stewart_Moce1_Approch_From_Outside_To_Inside
