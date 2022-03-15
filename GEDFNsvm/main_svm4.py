## Yunchuan Kong
## 2018 Copyright Reserved

## 安装 tensorflow, https://blog.csdn.net/u012270544/article/details/96424907


from __future__ import print_function
import numpy as np
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

from sklearn.utils import shuffle
from sklearn import metrics
from random import seed
import time
import sys


import pandas as pd


tf.reset_default_graph()
# =============================================================================
# file1 = sys.argv[1]
# file2 = sys.argv[2]
# out_file = sys.argv[3]
# =============================================================================

#import os
#os.chdir('D:/E/博士/DR_paper/Paper2/参考文献/新查阅的AIreview文献/2018_A graph-embedded deep/GEDFN-master')

# sample * gene
file_train = "TCGA_pro_outcome_TN_log_comp_allgeneT.txt"
file_test  = "TCGA_pro_outcome_TN_log_comp_allgeneT.txt"
file2 = "adjmatrix_comp_allgene.csv"


## load in data
partition = pd.read_table(file2, sep=',', index_col=0)    # 16144*16144 
expression_train = pd.read_table(file_train, sep='\t', index_col=0)    # 36*127
expression_test  = pd.read_table(file_test,  sep='\t', index_col=0)    # 12*127


label_train = np.array(expression_train['outcome'], dtype=int)    # 36
label_test  = np.array(expression_test['outcome'],  dtype=int)    # 12


labels_train = []
for l in label_train:
    if l == 1:
        labels_train.append([0,1])
    else:
        labels_train.append([1,0])
labels_train = np.array(labels_train,dtype=int)


labels_test = []
for l in label_test:
    if l == 1:
        labels_test.append([0,1])
    else:
        labels_test.append([1,0])
labels_test = np.array(labels_test,dtype=int)


exp_train = expression_train.drop('outcome', axis=1)  # 36*126
exp_test  = expression_test.drop('outcome',  axis=1)  # 12*126


x_train = exp_train  # 36*126
x_test  = exp_test  # 12*126
y_train = labels_train  # 36*2
y_test  = labels_test  # 123*2



## hyper-parameters and settings
L2 = False
max_pooling = False
droph1 = False
learning_rate = 0.0001
training_epochs = 100
batch_size = 8
display_step = 1

## the constant limit for feature selection
gamma_c = 50
gamma_numerator = np.sum(partition, axis=0)    # 按列求和
gamma_denominator = np.sum(partition, axis=0)

## Notice！！！   example is 

#gamma_numerator[np.where(gamma_numerator>gamma_c)] = gamma_c   
#gamma_numerator[gamma_numerator[gamma_numerator>gamma_c]] = gamma_c
np.where(gamma_numerator>gamma_c)
gamma_numerator[33]
gamma_numerator[33] = gamma_c



n_hidden_1 = np.shape(partition)[0]
n_hidden_2 = 64
n_hidden_3 = 16
n_classes = 2
n_features = np.shape(exp_train)[1]  # 100

## initiate training logs
loss_rec = np.zeros([training_epochs, 1])    # 1 列
training_eval = np.zeros([training_epochs, 2])    # 2 列

def max_pool(mat): ## input {mat}rix

    def max_pool_one(instance):
        return tf.reduce_max(tf.multiply(tf.matmul(tf.reshape(instance, [n_features, 1]), tf.ones([1, n_features]))
                                         , partition)
                             , axis=0)

    out = tf.map_fn(max_pool_one, mat, parallel_iterations=1000, swap_memory=True)
    return out

def multilayer_perceptron(x, weights, biases, keep_prob):
    layer_1 = tf.add(tf.matmul(x, tf.multiply(weights['h1'], partition)), biases['b1'])
    layer_1 = tf.nn.relu(layer_1)
    if max_pooling:
        layer_1 = max_pool(layer_1)
    if droph1:
        layer_1 = tf.nn.dropout(layer_1, keep_prob=keep_prob)

    layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
    layer_2 = tf.nn.relu(layer_2)
    layer_2 = tf.nn.dropout(layer_2, keep_prob=keep_prob)

    layer_3 = tf.add(tf.matmul(layer_2, weights['h3']), biases['b3'])
    ## Do not use batch-norm
    # layer_3 = tf.contrib.layers.batch_norm(layer_3, center=True, scale=True,
    #                                   is_training=is_training)
    layer_3 = tf.nn.relu(layer_3)
    layer_3 = tf.nn.dropout(layer_3, keep_prob=keep_prob)

    out_layer = tf.matmul(layer_3, weights['out']) + biases['out']
    return out_layer


x = tf.placeholder(tf.float32, [None, n_features])
y = tf.placeholder(tf.int32, [None, n_classes])
keep_prob = tf.placeholder(tf.float32)
lr = tf.placeholder(tf.float32)

weights = {
    'h1': tf.Variable(tf.truncated_normal(shape=[n_features, n_hidden_1], stddev=0.1)),
    'h2': tf.Variable(tf.truncated_normal(shape=[n_hidden_1, n_hidden_2], stddev=0.1)),
    'h3': tf.Variable(tf.truncated_normal(shape=[n_hidden_2, n_hidden_3], stddev=0.1)),
    'out': tf.Variable(tf.truncated_normal(shape=[n_hidden_3, n_classes], stddev=0.1))

}

biases = {
    'b1': tf.Variable(tf.zeros([n_hidden_1])),
    'b2': tf.Variable(tf.zeros([n_hidden_2])),
    'b3': tf.Variable(tf.zeros([n_hidden_3])),
    'out': tf.Variable(tf.zeros([n_classes]))
}

# Construct model
pred = multilayer_perceptron(x, weights, biases, keep_prob)

# Define loss and optimizer
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
if L2:
    reg = tf.nn.l2_loss(weights['h1']) + tf.nn.l2_loss(weights['h2']) + \
          tf.nn.l2_loss(weights['h3']) + tf.nn.l2_loss(weights['out'])
    cost = tf.reduce_mean(cost + 0.01 * reg)
optimizer = tf.train.AdamOptimizer(learning_rate=lr).minimize(cost)

## Evaluation
correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
y_score = tf.nn.softmax(logits=pred)

var_left = tf.reduce_sum(tf.abs(tf.multiply(weights['h1'], partition)), 0)
var_right = tf.reduce_sum(tf.abs(weights['h2']), 1)
var_importance = tf.add(tf.multiply(tf.multiply(var_left, gamma_numerator), 1./gamma_denominator), var_right)

with tf.Session() as sess:

    sess.run(tf.global_variables_initializer())
    total_batch = int(np.shape(x_train)[0] / batch_size)

    ## Training cycle
    for epoch in range(training_epochs):
        avg_cost = 0.
        x_tmp, y_tmp = shuffle(x_train, y_train)
        # Loop over all batches
        for i in range(total_batch-1):
            batch_x, batch_y = x_tmp[i*batch_size:i*batch_size+batch_size], \
                                y_tmp[i*batch_size:i*batch_size+batch_size]

            _, c= sess.run([optimizer, cost], feed_dict={x: batch_x, y: batch_y,
                                                        keep_prob: 0.9,
                                                        lr: learning_rate
                                                        })
            # Compute average loss
            avg_cost += c / total_batch

        del x_tmp
        del y_tmp

        ## Display logs per epoch step
        if epoch % display_step == 0:
            loss_rec[epoch] = avg_cost
            acc, y_s = sess.run([accuracy, y_score], feed_dict={x: x_train, y: y_train, keep_prob: 1})
            auc = metrics.roc_auc_score(y_train, y_s)
            training_eval[epoch] = [acc, auc]
            print ("Epoch:", '%d' % (epoch+1), "cost =", "{:.9f}".format(avg_cost),
                    "Training accuracy:", round(acc,3), " Training auc:", round(auc,3))

        if avg_cost <= 0.1:
            print("Early stopping.")
            break

    ## Testing cycle
    acc, y_s = sess.run([accuracy, y_score], feed_dict={x: x_test, y: y_test, keep_prob: 1})
    auc = metrics.roc_auc_score(y_test, y_s)
    var_imp = sess.run([var_importance])
    var_imp = np.reshape(var_imp, [n_features])
    print("*****=====", "Testing accuracy: ", acc, " Testing auc: ", auc, "=====*****")

np.savetxt("var_impo_svm4.csv", var_imp, delimiter=",")



















