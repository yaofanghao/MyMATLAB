'''
Author: yao fanghao
Date: 2023-12-11 21:15:11
LastEditTime: 2023-12-12 15:36:56
LastEditors: yao fanghao
'''
import numpy as np
import torch
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

df = pd.read_excel('all_label.xlsx')

train_data = df.iloc[:, :8].values  # 前八列
train_labels = df.iloc[:, 8].values  # 第九列

train_data = torch.tensor(train_data, dtype=torch.float32)
train_labels = torch.tensor(train_labels, dtype=torch.long)

class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(NeuralNetwork, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out

input_size = 8  # 输入特征
hidden_size = 128  # 隐藏层的大小
num_classes = 9  # 输出类别

learning_rate = 0.001  # 学习率
num_epochs = 20  # 训练轮数

model = NeuralNetwork(input_size, hidden_size, num_classes)
criterion = nn.CrossEntropyLoss()
optimizer = optim.SGD(model.parameters(), lr=learning_rate)

total_steps = len(train_data)
for epoch in range(num_epochs):
    for i in range(total_steps):
        x = train_data[i].unsqueeze(0) 
        labels = train_labels[i].unsqueeze(0)

        outputs = model(x)
        loss = criterion(outputs, labels)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch+1}/{num_epochs}], Step [{i+1}/{total_steps}], Loss: {loss.item()}')

data_test = pd.read_excel("test_label.xlsx")
test_x = np.array(data_test.iloc[:, 0:8])
test_y = np.array(data_test.iloc[:, 8])
x_test = torch.FloatTensor(test_x)
y_test = torch.LongTensor(test_y)

# 进行准确率评估
output_prediction = model(x_test)
output_pred = torch.argmax(output_prediction, dim=1)
output_pred = output_pred.numpy()
acc_test = 0
for item in range(len(test_y)):
    if output_pred[item] == test_y[item]:
        acc_test += 1
acc_test /= len(test_y)
acc_test = round(acc_test, 3)
print('test of acc:{}'.format(acc_test))