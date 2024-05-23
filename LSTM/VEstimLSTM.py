# Author: Junran Chen
# Date: 2023-June-23
# Function: module for LSTM model
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
from datetime import datetime, timedelta


class VEstimLSTM(nn.Module):
    def __init__(self, hidden_size, input_size, layers):
        super().__init__()
        self.lstm = nn.LSTM(input_size=input_size, hidden_size=hidden_size,
                            num_layers=layers, batch_first=True)
        self.linear = nn.Linear(hidden_size, 1)
        # self.h_s = None
        # self.h_c = None

    def forward(self, x, h_s, h_c):
        # The h_s, h_c is defaulted to 0 every time, so only remember last 500-second data
        y, (h_s, h_c) = self.lstm(x, (h_s, h_c))
        y = self.linear(y)
        return y, (h_s, h_c)