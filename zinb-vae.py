######
#  ZINB-VAE Code
#  Romain Lopez, UC Berkeley, 2017
#  TF 1.3.0
#  python 2.7.12
#  
#  this code is not meant to work as such 
#  but to show the computational graph and optimization process
#######

import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import numpy as np
import tensorflow as tf
from tensorflow.contrib import slim
from tensorflow.contrib.framework import arg_scope
import time
import pickle
import h5py
from ZIAE_helper import *
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

####
# Helpful Functions
####
# Gaussian Sample Layer in TF
def gaussian_sample(mean, var, scope=None):
    with tf.variable_scope(scope, 'gaussian_sample'):
        sample = tf.random_normal(tf.shape(mean), mean, tf.sqrt(var))
        sample.set_shape(mean.get_shape())
        return sample

# Dense Layer in TF Slim
def dense(x,
          num_outputs,
          STD=0.01,
          lambda_L2=0.001,
          keep_prob=None,
          scope=None,
          activation=None,
          reuse=None,
          bn=False,
          phase=None):
    
    if keep_prob is not None:
        output = tf.layers.dropout(x, rate=keep_prob, training=phase)
    else:
        output = tf.identity(x)
        
    output = slim.fully_connected(output, num_outputs, activation_fn=None, \
                                weights_initializer=tf.truncated_normal_initializer(stddev=STD), \
                                  weights_regularizer=slim.l2_regularizer(lambda_L2)
                                 )
    if bn:
        output = tf.layers.batch_normalization(output, training=phase)
    if activation: 
        output = activation(output)
    return output

def log_zinb_positive(x, mu, theta, pi, zi=True, eps=1e-8):
    """
    log likelihood of a batch according to a zinb model
    =====
    mu: mean (>0)
    theta: inverse dispersion parameter (>0)
    pi:logit of the dropout parameter 
    """
    case_zero = tf.nn.softplus(- pi + theta * tf.log(theta + eps) - theta * tf.log(theta + mu + eps)) \
                                - tf.nn.softplus( - pi)
    case_non_zero = - pi - tf.nn.softplus(- pi) \
                                + theta * tf.log(theta + eps) - theta * tf.log(theta + mu + eps) \
                                + x * tf.log(mu + eps) - x * tf.log(theta + mu + eps) \
                                + tf.lgamma(x + theta) - tf.lgamma(theta) - tf.lgamma(x + 1)
    mask = tf.cast(tf.less(x, eps), tf.float32)
    res = tf.multiply(mask, case_zero) + tf.multiply(1 - mask, case_non_zero)
    return tf.reduce_sum(res, axis=-1)

####
# Setting of the script
####

#Constant parameters, not to be changed
STD = 0.01 # neural nets weight initialization
n_input = 720 # number of genes
M = 128  # batch size during training
n_latent = 10 # dimension of the latent space

# Parameters to be changed by re-running the whole script
n_layers = 2
GENERATIVE_MODEL = True # learn a AE or a VAE
num_epochs = 300 # from 10 for the 1M to 150 for the small dataset
n_hidden = 128 # dimension for hidden layers: usually 128 or 256
complex_dispersion = True # full NN if true, one parameter per gene otherwise
stochastic_batches = False # Use False for a hdf5 file

# Dynamic parameters
dropout_nns = 0.1 # 0.1 to 0.5
learning_rates = 0.005 # 0.005 to 0.0001
epsilons = 0.01 # 1e-1  to 1e-4

# Import your own dataset
dataset = "huge_1M"
sample_size = 10000
X_train, X_test, X_train_log, X_test_log, c_train, c_test \
        = load_dataset(dataset,  prefix="/home/ubuntu/data", \
                       n_cell=sample_size)

# load test set in RAM
X_test = X_test[:]
print(X_train.shape, "learning_set")
print("DATA LOADED")


####
# Now this is the definition of the model in TF
####


# placeholders
tf.reset_default_graph()
x = tf.placeholder(tf.float32, (None, n_input), name='x')
training_phase = tf.placeholder(tf.bool, (), name='training_phase')
keep_prob = tf.placeholder(tf.float32, ())
l2_weight = tf.placeholder(tf.float32, ())
    
# variational distribution or encoder
with tf.variable_scope('variational_gaussian'):
    # q_z | x
    h = dense(x, n_hidden, scope='layer1', activation=tf.nn.relu, \
                STD=STD, bn=True, keep_prob=keep_prob, phase=training_phase)
    for layer in range(2, n_layers + 1):
        h = dense(h, n_hidden, scope='layer' + str(layer), activation=tf.nn.relu, \
                STD=STD, bn=True, keep_prob=keep_prob, phase=training_phase)

    if GENERATIVE_MODEL:
        qz_m = dense(h, n_latent, scope='mean', activation=None, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
        qz_v = dense(h, n_latent, scope='var', activation=tf.exp, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
    else:
        z  = dense(h, n_latent, scope='z', activation=None, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
        
if GENERATIVE_MODEL:
    with tf.variable_scope('sample'):
        z = gaussian_sample(qz_m, qz_v)
    
    
    
# generative model or decoder   
with tf.variable_scope('generative_model'):
    # p_x | z
    h = dense(z, n_hidden, scope='layer1', activation=tf.nn.relu, \
                STD=STD, bn=True, keep_prob=keep_prob, phase=training_phase)
    for layer in range(2, n_layers + 1):
        h = dense(h, n_hidden, scope='layer' + str(layer), activation=tf.nn.relu, \
                STD=STD, bn=True, keep_prob=keep_prob, phase=training_phase)
    
    px_rate = dense(h, n_input, scope="px_rate", activation=tf.exp, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
    
    px_dropout = dense(h, n_input, scope='px_dropout', activation=None, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
    if complex_dispersion:
        px_r = dense(h, n_input, scope='px_r', activation=None, \
                STD=STD, bn=False, keep_prob=None, phase=training_phase)
    else:
        px_r = tf.Variable(tf.random_normal([n_input]), name="r")

# variational lower bound or reconstruction error
with tf.name_scope('loss'):
    recon = log_zinb_positive(x, px_rate, tf.exp(px_r), px_dropout, zi=True)
    if GENERATIVE_MODEL:
        kl_gauss = 0.5 * tf.reduce_sum(tf.square(qz_m) + qz_v - tf.log(1e-8 + qz_v) - 1, 1)
        ELBO_gau = tf.reduce_mean(recon - kl_gauss)
        loss_gau = - ELBO_gau
    else:
        loss_gau = - tf.reduce_mean(recon)

# CLUSTERING MODEL: See Variational Deep Embedding Paper (IJCAI2017)
K = 7
if GENERATIVE_MODEL:
    with tf.variable_scope('clustering_model'):
        pi = tf.Variable(tf.ones([K])/K, trainable=False)
        mu = tf.Variable(tf.random_normal([K, n_latent]), name="mu")
        var = tf.Variable(tf.nn.softplus(tf.random_normal([K, n_latent])), name="sigma", trainable=True)

    with tf.variable_scope('variational_cat'):
        # q c | x = p(c | z)
        z_repeat = tf.reshape(tf.tile(z, [1, K]), [-1, K, n_latent])
        qz_m_repeat = tf.reshape(tf.tile(qz_m, [1, K]), [-1, K, n_latent])
        qz_v_repeat = tf.reshape(tf.tile(qz_v, [1, K]), [-1, K, n_latent])

        pc_logits = tf.log(pi + 1e-8) + tf.reduce_sum(- 0.5 * tf.square(z_repeat - mu) / var \
                                        - 0.5 * tf.log(var + 1e-8) - 0.5 * tf.log(2 * np.pi), 2) + 1e-5
        pc_logits = tf.identity(pc_logits, "pc_logits")
        pc = tf.nn.softmax(pc_logits, name="pc")

    with tf.name_scope('clustering_loss'):

        kl_gaussian_mixture = 0.5 * tf.reduce_sum(pc * tf.reduce_sum(tf.square(qz_m_repeat - mu) / var \
                        + qz_v_repeat / var \
                        + tf.log(1e-8 + var) - tf.log(1e-8 + qz_v_repeat) - 1, 2), 1)
        kl_cat = tf.reduce_sum(pc * ( tf.log(pc + 1e-8) - tf.log(pi + 1e-8)), 1)
        ELBO_gmm = tf.reduce_mean(recon - kl_gaussian_mixture - kl_cat)
        loss_gmm = - ELBO_gmm



print("TF MODEL LOADED")

# Complete Optimizer to learn the model
update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, epsilon=epsilon)
with tf.control_dependencies(update_ops):
    train_step = optimizer.minimize(loss_gau)
    
# Test time optimizer in the AE case to compute validation log likelihood
update_ops_test = tf.get_collection(tf.GraphKeys.UPDATE_OPS, "variational_gaussian")
test_vars = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "variational_gaussian")
optimizer_test = tf.train.AdamOptimizer(learning_rate=0.001, epsilon=0.1)
with tf.control_dependencies(update_ops_test):
    test_step = optimizer_test.minimize(loss_gau, var_list=test_vars)

if GENERATIVE_MODEL:
    # GMM prior optimizer in the labels application
    optimizer_gmm = tf.train.AdamOptimizer(learning_rate=learning_rate, epsilon=epsilon)
    with tf.control_dependencies(update_ops):
        gmm_step = optimizer_gmm.minimize(loss_gmm)

# Session creation
sess = tf.Session()
sess.run(tf.global_variables_initializer())

####
# TRAINING PROCEDURE
####

for t in range(iterep * num_epochs):
    
    # arange data in batches
    if stochastic_batches:
        x_train = X_train[np.random.choice(np.arange(X_train.shape[0]), size=M)].astype(np.float32)
    else:
        x_train = X_train[M * (t%iterep) :M * (t%iterep + 1) ].astype(np.float32) 
    x_test = X_test[np.random.choice(np.arange(X_test.shape[0]), size=M)].astype(np.float32)

    #prepare data dictionaries    
    dic_train = {x: x_train, training_phase:True, keep_prob:dropout_nn}
    dic_test = {x: x_test, training_phase:False, keep_prob:dropout_nn}
    
    # run an optimization step
    _, l = sess.run([train_step, loss_gau], feed_dict=dic_train)
    end_epoch, epoch = t % iterep == 0, t / iterep
    
    if end_epoch:
        now = time.time()
        print epoch
        
        l_tr = sess.run((loss_gau), feed_dict=dic_train)
        l_t = sess.run((loss_gau), feed_dict=dic_test)
        print 'Train / Test performance:', l_tr, l_t

dic_full = {x: X_test, training_phase:False, keep_prob:dropout_nn}

#if there are some labels
#z_test = sess.run(z, feed_dict=dic_full)
#sil_score = silhouette_score(z_test, c_test)
#print sil_score

# if interested in a likelihood
l = sess.run((loss_gau), feed_dict=dic_full)
print("likelihood on the test set:", l)

####
#   Optimization procedure of the encoder on the test set
####

for t in range(iterep * num_epochs):
    

    x_test = X_test[np.random.choice(np.arange(X_test.shape[0]), size=100)].astype(np.float32)

    #prepare data dictionaries
    dic_test = {x: x_test, training_phase:False, keep_prob:dropout_nn}
    dic_full = {x: X_test, training_phase:False, keep_prob:dropout_nn} 
    
    # run an optimization set
    _, l = sess.run([test_step, loss_gau], feed_dict=dic_test)
    end_epoch, epoch = t % iterep == 0, t / iterep
    
    if end_epoch:
        now = time.time()
        print epoch
        
        l_tr = sess.run((loss_gau), feed_dict=dic_test)
        l_t = sess.run((loss_gau), feed_dict=dic_full)
        print 'Train / Test performance:', l_tr, l_t

####
# Optimization procedure for the GMM PRIOR
####

# project on the latent space and initialize a GMM to send the parameters to the GMM PRIOR
z_test = sess.run(qz_m, feed_dict={'x:0': X_test, training_phase:False, keep_prob:dropout_nn})
z_ = sess.run(qz_m, feed_dict={'x:0': X_train, training_phase:False, keep_prob:dropout_nn})
model = GaussianMixture(n_components=K, covariance_type="diag")
model.fit(z_)
init_res, _ = cluster_acc(model.predict(z_test), c_test)
sess.run(mu.assign(model.means_))
sess.run(var.assign(model.covariances_))

iterep = 500
for i in range(iterep * 100):
        x_train = X_train[np.random.choice(np.arange(X_train.shape[0]), size=M)].astype(np.float32)
        x_test = X_test[np.random.choice(np.arange(X_test.shape[0]), size=M)].astype(np.float32)
        train_dic = {x: x_train, training_phase:True, keep_prob:dropout_nn}
        test_dic = {x: x_test, training_phase:False, keep_prob:dropout_nn}
        _, l = sess.run([gmm_step, loss_gmm], feed_dict=train_dic)
        
        end_epoch, epoch = i % iterep == 0, i / iterep
        if end_epoch:
            l_tr = sess.run((loss_gmm), feed_dict=train_dic)
            l_t, z_test = sess.run((loss_gmm, qz_m), feed_dict=test_dic)
            print 'Train / Test performance:', l_tr, l_t

full_dic = {x: X_test, training_phase:False, keep_prob:dropout_nn}
z_test = sess.run(qz_m, feed_dict=full_dic)
print "Silhouette TEST", silhouette_score(z_test, c_test)
