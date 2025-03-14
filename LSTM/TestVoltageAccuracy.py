# Author: Junran Chen
# Date: 2023-June-23
# Function: Test the voltage estimation accuracy.
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import timeit
import torch.optim as optim
import torch.utils.data as data
from datetime import datetime, timedelta
from VEstimLSTM import VEstimLSTM
# ---------------------Tested on sequence input-------------------------------------------
# BATCH_SIZE = 3
# HIDDEN_SIZE = 16
# loss_fn = nn.MSELoss()
# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')  # select CPU or GPU
# [X_Test, Y_Test] = torch.load('./data/Combined_Testing15_2sequnce-May-2023.pt')
# testLoader = data.DataLoader(data.TensorDataset(X_Test, Y_Test), shuffle=False,
#                              batch_size=BATCH_SIZE, drop_last=True)
# # best_model = VEstimLSTM(hidden_size=HIDDEN_SIZE)
# best_model = torch.load('./result/LSTM_June24.model', map_location='cuda:0').get('model').cuda()
# best_model.eval()
# y_predict = []
# y_measured = []
# total_valid_loss = []
# with torch.no_grad():
#     for X_batch, y_batch in testLoader:
#         h_s = torch.ones(2, BATCH_SIZE, 16)  # 2-layers, 100-batch, 16-hidden layers.
#         h_c = torch.zeros(2, BATCH_SIZE, 16)
#         h_s = h_s.to(device)
#         h_c = h_c.to(device)
#         X_batch = X_batch.to(device)
#         y_batch = y_batch.to(device)
#         y_pred, (h_s, h_c) = best_model(X_batch, h_s, h_c)
#         test_rmse = torch.sqrt(loss_fn(y_pred[:, -1, :], y_batch))
#         total_valid_loss.append(test_rmse.item())
#         y_measured = y_measured + y_batch.cpu().numpy().tolist()
#         y_predict = y_predict + y_pred[:, -1, :].cpu().numpy().tolist()
# print('MSE=', np.mean(total_valid_loss))
# y_measured = np.asarray(y_measured)
# y_predict = np.asarray(y_predict)
# compareV = np.concatenate((y_predict, y_measured), axis=1)
# np.savetxt("./result/SOP_25degC_June_24.csv", compareV)

# ---------------------------Tested on data without sequence---------------------------------
# loss_fn = nn.MSELoss()
# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')  # select CPU or GPU
# [X_Test, Y_Test] = torch.load('./data/_25degC_US06.pt')
# X_Test = torch.tensor(X_Test).cuda()
# best_model = torch.load('./result/LSTM_Aug18_CUDA1.model', map_location='cuda:0').get('model').cuda()
# best_model.eval()
# y_predict = []
# y_measured = []
# total_valid_loss = []
# h_s = torch.zeros(2, 16)  # 2-layers, 100-batch, 16-hidden layers.
# h_c = torch.zeros(2, 16)
# with torch.no_grad():
#     h_s = h_s.to(device)
#     h_c = h_c.to(device)
#     y_pred, (h_s, h_c) = best_model(X_Test, h_s, h_c)
#     y_pred = np.asarray(y_pred.cpu())
# y_pred = np.array(y_pred)
# Y_Test = np.array(Y_Test)
# error_P = (y_pred - Y_Test)
# RMSEP = np.sqrt(np.mean(np.square(error_P)))*1000
# print(RMSEP)
# compareV = np.concatenate((y_pred, Y_Test), axis=1)
# np.savetxt("./result/Voltage_test.csv", compareV)
# ---------------------------Tested on data from csv file---------------------------------

loss_fn = nn.MSELoss()
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')  # select CPU or GPU
test = pd.read_csv('./Samsung30T_driveCycle.csv')
Y_Test = test[["V"]].values.astype('float32')
test = test[["Ah", "T", "P"]].values.astype('float32')
test[:, 0] = (   # align data to fully charge
                            test[:, 0] + 3 - np.max(test[:, 0]))
X_Test = torch.tensor(test[26802:29693]).cuda()
Y_Test = torch.tensor(Y_Test[26802:29693])
best_model = torch.load('./LSTM_Aug18_CUDA1.model', map_location='cuda:0').get('model').cuda()
best_model.eval()
y_predict = []
y_measured = []
total_valid_loss = []
start = timeit.default_timer()
h_s = torch.zeros(2, 16)  # 2-layers, 100-batch, 16-hidden layers.
h_c = torch.zeros(2, 16)
with torch.no_grad():
    h_s = h_s.to(device)
    h_c = h_c.to(device)
    y_pred, (h_s, h_c) = best_model(X_Test, h_s, h_c)
    # i = 1
    # while i < len(X_Test):
    #     y_pred, (h_s, h_c) = best_model(torch.unsqueeze(X_Test[i, :], 0), h_s, h_c)
    #     i = i+1
    y_pred = np.asarray(y_pred.cpu())
y_pred = np.array(y_pred) - 0.15
Y_Test = np.array(Y_Test)
error_P = (y_pred[999:-1] - Y_Test[999:-1])
RMSEP = np.sqrt(np.mean(np.square(error_P)))*1000
print(RMSEP)
stop = timeit.default_timer()
print('Time: ', stop - start)
compareV = np.concatenate((y_pred, Y_Test), axis=1)
np.savetxt("./Voltage_test.csv", compareV)

