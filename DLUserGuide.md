# 部分学习路线
# 参考资料：官方deeping learning toolbox **用户指南**
用户指南.pdf 网盘链接： 

链接：https://pan.baidu.com/s/1-w0DU2dDJ-wmNokcPWgZnA 提取码：1234

数字表示在用户指南中的页码

# 变量管理
* save mydata a x
* load mydata

# **深度学习层列表** 
非常详细的介绍 1-18--1-23
输入，卷积和全连接，序列，激活，归一化、丢弃和裁剪，池化和去池化，组合，目标检测，生成对抗，输出层

# 导入与训练网络
从Tensorflow-Keras 1-13 
* 可以导入.h5文件格式的网络和权重
* 参阅importKerasNetwork
  
在ONNX中导入和导出网络 1-13
* 以ONNX为中间格式，在Tensorflow、MATLAB等框架之间交互模型

# CNN
非常详细的原理介绍 1-25--1-32

设置参数并训练 1-35--1-37
* 保存检查点：使用trainingOptions中的CheckpointPath，双击或命令行load加载文件；
* 有关示例参考Resume Training from Checkpoint Network

# deep network designer工具箱
迁移学习 2-2--2-14

构建网络 2-15--2-20

-------------------------------

# 图像深度学习
详细介绍：
https://ww2.mathworks.cn/help/deeplearning/deep-learning-with-images.html

使用摄像头 3-2--3-5

训练网络 3-6--3-12

使用残差网络 3-13--3-21

使用GoogLeNet 3-22--3-26

使用预训练网络提取特征 3-27--3-31

使用AlexNet迁移学习 3-32--3-38

创建简单网络 3-39--3-43

拟合回归模型预测手写数字旋转角度 3-44--3-51

将分类转化为回归网络 3-52--3-56

生成对抗网络GAN 3-57--3-67

训练变分自编码器VAE 3-68--3-78

------------------------------

# 深度学习调整和可视化
详细介绍：
https://ww2.mathworks.cn/help/deeplearning/deep-learning-tuning-and-visualization.html?s_tid=CRUX_lftnav

GoogLeNet的Deep Dream 5-2--5-7

贝叶斯优化 5-8--5-16

并行训练 5-17--5-21

监控训练进度 5-22--5-25

自定义训练期间的输出 5-26--5-29

热力图CAM可视化预测过程 5-30--5-34

可视化CNN激活特征 5-35--5-45

可视化CNN特征 5-46--5-52

------------------------------

#计算机视觉
Augment Bounding Boxes for Object Detection
* https://ww2.mathworks.cn/help/deeplearning/ug/bounding-box-augmentation-using-computer-vision-toolbox.html

YOLOv2 8-2--8-10

DeeplabV3+的语义分割 8-11--8-25

扩张卷积的语义分割 8-26--8-30

UNet的多光谱语义分割 8-31--8-47

定义使用Tversky Loss的像素分类层 8-59--8-65

RCNN的目标检测 8-66--8-77

Faster RCNN的目标检测 8-78--8-86

Grad-CAM可视化
* https://ww2.mathworks.cn/help/deeplearning/ug/explore-semantic-segmentation-network-using-gradcam.html

------------------------------

# 图像处理
图像增强
* https://ww2.mathworks.cn/help/deeplearning/ug/image-augmentation-using-image-processing-toolbox.html

使用预训练网络去除图像噪声 9-2--9-7

单图像超分辨率 9-8--9-20

DnCNN的图像去块 9-21--9-32

CAN的图像处理算子逼近 9-33--9-46

VGG19作为预训练的风格迁移 9-47--9-55

------------------------------

#自动驾驶

Faster RCNN的车辆检测 10-2--10-11

使用单目相机和语义分割创建占据栅格 10-12--10-25

语义分割自动标注Ground Truth标签
* https://ww2.mathworks.cn/help/driving/ug/automate-ground-truth-labeling-for-semantic-segmentation.html

------------------------------

# 导入、导出和自定义
从外部深度学习平台导入和导出网络
* https://ww2.mathworks.cn/help/deeplearning/deep-learning-import-and-export.html

自定义深度学习层 17-2--17-15
* 层模板、中间层、检查层有效性、在网络中包含层、输出层

自定义训练循环、损失函数和网络 17-33--37
* 训练循环、损失函数、自动微分
* https://ww2.mathworks.cn/help/deeplearning/deep-learning-custom-training-loops.html

权重初始化 17-16--17-27
* 定义和比较

基于预训练的Keras层组合网络 17-28--17-32

--------------------------------

# 数据预处理











------------------------------


