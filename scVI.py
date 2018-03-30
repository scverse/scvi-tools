"""
Code for the paper single-cell Variational Inference (scVI) paper

Romain Lopez, Jeffrey Regier, Michael Cole, Michael Jordan, Nir Yosef
EECS, UC Berkeley

"""


import functools
import tensorflow as tf
import numpy as np
from tensorflow.contrib import slim


def dense(x, 
          num_outputs,
          STD=0.01,
          keep_prob=None,
          activation=None,
          bn=False,
          phase=None):
    """
    Defining the elementary layer of the network

    Note:
    We adjust the standard deviation of the weight initialization to 0.01. This is useful considering the large counts in the data and guarantees numerical stability of the algorithm.  
    Batchnorm paper: https://arxiv.org/abs/1502.03167
    Dropout paper: http://jmlr.org/papers/v15/srivastava14a.html

    Variables:
    x: tensorflow variable
    num_outputs: number of outputs neurons after the dense layer
    keep_prob: float number for probability of keeping an individual neuron for dropout layer
    activation: tensorflow activation function (tf.exp, tf.nn.relu...) for this layer
    bn: bool to use batchnorm for this layer
    phase: tensorflow boolean node indicating whether training of testing phase (see dropout and batchnorm paper)
    """
    output = tf.identity(x)
    
    if keep_prob is not None:
        output = tf.layers.dropout(output, rate=keep_prob, training=phase)

    output = slim.fully_connected(output, num_outputs, activation_fn=None, \
                            weights_initializer=tf.truncated_normal_initializer(stddev=STD))

    if bn:
        output = tf.layers.batch_normalization(output, training=phase)
        
    if activation: 
        output = activation(output)
            
    return output



def gaussian_sample(mean, var, scope=None):
    """
    Function to sample from a multivariate gaussian with diagonal covariance in tensorflow

    Note:
    This layer can either be parametrized by the variance or the log variance in a variational autoencoder. 
    We found by trials that it does not matter much

    Variables:
    mean: tf variable indicating the minibatch mean (shape minibatch_size x latent_space_dim)
    var: tf variable indicating the minibatch variance (same shape)
    """
    with tf.variable_scope(scope, 'gaussian_sample'):
        sample = tf.random_normal(tf.shape(mean), mean, tf.sqrt(var))
        sample.set_shape(mean.get_shape())
        return sample

    

def log_zinb_positive(x, mu, theta, pi, eps=1e-8):
    """
    log likelihood (scalar) of a minibatch according to a zinb model. 

    Notes:
    We parametrize the bernouilli using the logits, hence the softplus functions appearing
    
    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    pi: logit of the dropout parameter (real support) (shape: minibatch x genes)
    eps: numerical stability constant
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

def log_nb_positive(x, mu, theta, eps=1e-8):
    """
    log likelihood (scalar) of a minibatch according to a nb model. 
    
    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    eps: numerical stability constant
    """    
    res = tf.lgamma(x + theta) - tf.lgamma(theta) - tf.lgamma(x + 1) + x * tf.log(mu + eps) \
                                - x * tf.log(theta + mu + eps) + theta * tf.log(theta + eps) \
                                - theta * tf.log(theta + mu + eps)
    return tf.reduce_sum(res, axis=-1)
    

def doublewrap(function):
    """
    A decorator decorator, allowing to use the decorator to be used without
    parentheses if not arguments are provided. All arguments must be optional.

    Notes:
    https://gist.github.com/danijar/8663d3bbfd586bffecf6a0094cd116f2
    """
    @functools.wraps(function)
    def decorator(*args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            return function(args[0])
        else:
            return lambda wrapee: function(wrapee, *args, **kwargs)
    return decorator

@doublewrap
def define_scope(function, scope=None, *args, **kwargs):
    """
    A decorator for functions that define TensorFlow operations. The wrapped
    function will only be executed once. Subsequent calls to it will directly
    return the result so that operations are added to the graph only once.
    The operations added by the function live within a tf.variable_scope(). If
    this decorator is used with arguments, they will be forwarded to the
    variable scope. The scope name defaults to the name of the wrapped
    function.

    Notes:
    https://gist.github.com/danijar/8663d3bbfd586bffecf6a0094cd116f2
    """
    attribute = '_cache_' + function.__name__
    name = scope or function.__name__
    @property
    @functools.wraps(function)
    def decorator(self):
        if not hasattr(self, attribute):
            with tf.variable_scope(name, *args, **kwargs):
                setattr(self, attribute, function(self))
        return getattr(self, attribute)
    return decorator

def mmd_fourier(x1, x2, bandwidth=2., dim_r=500):
    """
    Approximate RBF kernel by random features

    Notes:
    Reimplementation in tensorflow of the Variational Fair Autoencoder
    https://arxiv.org/abs/1511.00830
    """
    d = x1.get_shape().as_list()[1]
    rW_n = tf.sqrt(2. / bandwidth) * tf.random_normal([d, dim_r]) / np.sqrt(d)
    rb_u = 2 * np.pi * tf.random_uniform([dim_r])
    rf0 = tf.sqrt(2. / dim_r) * tf.cos(tf.matmul(x1, rW_n) + rb_u)
    rf1 = tf.sqrt(2. / dim_r) * tf.cos(tf.matmul(x2, rW_n) + rb_u)
    result = tf.reduce_sum((tf.reduce_mean(rf0, axis=0) - tf.reduce_mean(rf1, axis=0))**2)
    return tf.sqrt(result)

def mmd_rbf(x1, x2, bandwidths=1. / (2 * (np.array([1., 2., 5., 8., 10])**2))):
    """
    Return the mmd score between a pair of observations

    Notes:
    Reimplementation in tensorflow of the Variational Fair Autoencoder
    https://arxiv.org/abs/1511.00830
    """
    d1 = x1.get_shape().as_list()[1]
    d2 = x2.get_shape().as_list()[1]
    
    def K(x1, x2, gamma=1.): 
        dist_table = tf.expand_dims(x1, 0) - tf.expand_dims(x2, 1)
        return tf.transpose(tf.exp(-gamma * tf.reduce_sum(dist_table **2, axis=2)))

    # possibly mixture of kernels
    x1x1, x1x2, x2x2 = 0, 0, 0
    for bandwidth in bandwidths:
        x1x1 += K(x1, x1, gamma=np.sqrt(d1) * bandwidth) / len(bandwidths)
        x2x2 += K(x2, x2, gamma=np.sqrt(d2) * bandwidth) / len(bandwidths)
        x1x2 += K(x1, x2, gamma=np.sqrt(d1) * bandwidth) / len(bandwidths)

    return tf.sqrt(tf.reduce_mean(x1x1) - 2 * tf.reduce_mean(x1x2) + tf.reduce_mean(x2x2))

def mmd_objective(z, s, sdim):
    """
    Compute the MMD from latent space and nuisance_id

    Notes:
    Reimplementation in tensorflow of the Variational Fair Autoencoder
    https://arxiv.org/abs/1511.00830
    """
    
    #mmd_method = mmd_rbf
    mmd_method = mmd_fourier
    
    z_dim = z.get_shape().as_list()[1]

    # STEP 1: construct lists of samples in their proper batches
    z_part = tf.dynamic_partition(z, s, sdim)

                
    # STEP 2: add noise to all of them and get the mmd
    mmd = 0
    for j, z_j in enumerate(z_part):
        z0_ = z_j
        aux_z0 = tf.random_normal([1, z_dim])  # if an S category does not have any samples
        z0 = tf.concat([z0_, aux_z0], 0)
        if len(z_part) == 2:
            z1_ = z_part[j + 1]
            aux_z1 = tf.random_normal((1, z_dim))
            z1 = tf.concat([z1_, aux_z1], axis=0)
            return mmd_method(z0, z1)
        z1 = z
        mmd += mmd_method(z0, z1)
    return mmd


class scVIModel:

    def __init__(self, expression=None, batch_ind=None, num_batches=None, kl_scale=None, mmd_scale=None, phase=None,\
                 library_size_mean = None, library_size_var = None, apply_mmd=False, \
                 dispersion="gene", n_layers=1, n_hidden=128, n_latent=10, \
                 dropout_rate=0.1, log_variational=True, optimize_algo=None, zi=True):
        """
        Main parametrization of the scVI algorithm.

        Notes and disclaimer:
        + We recommend to put kl_scale to 1 for every tasks except clustering where 0 will lead better discrepency between the clusters
        + Applying a too harsh penalty will ruin your biology info. We recommend using less than a 100. From ongoing tests, using zero actually removes batch effects as well as the paper results.
        + We recommend the dispersion parameter to be gene specific (or batch-specific as in the paper) as in ZINB-WaVE if you do not have enough cells
        + To better remove library size effects between clusters, mention the log-library size prior for each batch (like in the paper)


        Variables:
        expression: tensorflow variable of shape (minibatch_size x genes), placeholder for input counts data
        batch_ind: tensorflow variable for batch indices (minibatch_size) with integers from 0 to n_batches - 1
        kl_scale: tensorflow variable for scalar multiplier of the z kl divergence
        mmd_scale: tensorflow variable for scalar multiplier of the MMD penalty
        phase: tensorflow variable for training phase
        library_size_mean = either a number or a list for each batch of the mean log library size
        library_size_var = either a number or a list for each batch of the variance of the log library size
        apply_mmd: boolean to choose whether to use a MMD penalty
        dispersion: "gene" (n_genes params) or "gene-batch" (n_genes x n_batches params) or "gene-cell" (a neural nets)
        n_layers: a integer for the number of layers in each neural net. We use 1 throughout the paper except on the 1M dataset where we tried (1, 2, 3) hidden layers
        n_hidden: number of neurons for each hidden layer. Always 128.
        n_latent: number of desired dimension for the latent space
        dropout_rate: rate to use for the dropout layer (see elementary layer function). always 0.1
        log_variational: whether to apply a logarithmic layer at the input of the variational network (for < 4000 cells datasets)
        optimize_algo: a tensorflow optimizer
        zi: whether to use a ZINB or a NB distribution
        """
        
        # Gene expression placeholder
        if expression is None:
            raise ValueError("provide a tensor for expression data")
        self.expression = expression
        
        print("Running scVI on "+ str(self.expression.get_shape().as_list()[1]) + " genes")
        self.log_variational = log_variational

        # batch correction
        if batch_ind is None:
            print("scVI will run without batch correction")
            self.batch = None
            self.apply_mmd = False
            
        else:
            if num_batches is None:
                raise ValueError("provide a comprehensive list of unique batch ids")
            self.batch_ind = batch_ind
            self.num_batches = num_batches
            self.batch = tf.one_hot(batch_ind, num_batches)
            self.mmd_scale = mmd_scale
            self.apply_mmd = apply_mmd

            print("Got " + str(num_batches) + "batches in the data")
            if self.apply_mmd:
                print("Will apply a MMD penalty")
            else: 
                print("Will not apply a MMD penalty")
        
        #kl divergence scalar
        if kl_scale is None:
            raise ValueError("provide a tensor for kl scalar")
        self.kl_scale = kl_scale
                
        #prior placeholder
        if library_size_mean is None or library_size_var is None:
            raise ValueError("provide prior for library size")
            
        if type(library_size_mean) in [float, np.float64] :
            self.library_mode = "numeric"
            self.library_size_mean = tf.to_float(tf.constant(library_size_mean))
            self.library_size_var = tf.to_float(tf.constant(library_size_var))
            
        else:
            if library_size_mean.get_shape().as_list()[0] != num_batches:
                raise ValueError("provide correct prior for library size (check batch shape)")
            else:
                self.library_mode = "list"
                self.library_size_mean = library_size_mean
                self.library_size_var = library_size_var 
                
        print("Will work on mode " + self.library_mode + " for incorporating library size")
        
        # high level model parameters
        if dispersion not in ["gene", "gene-batch", "gene-cell"]:
            raise ValueError("dispersion should be in gene / gene-batch / gene-cell")
        self.dispersion = dispersion
        
        print("Will work on mode " + self.dispersion + " for modeling inverse dispersion param")
        
        self.zi = zi
        if zi:
            print("Will apply zero inflation")
        
        # neural nets architecture
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.n_layers = n_layers
        self.n_input = self.expression.get_shape().as_list()[1] 

        print(str(self.n_layers) + " hidden layers at " + str(self.n_hidden) + " each for a final " + str(self.n_latent) + " latent space")
        
        # on training variables
        self.dropout_rate = dropout_rate
        if phase is None:
            raise ValueError("provide an optimization metadata (phase)")
        self.training_phase = phase
        if optimize_algo is None:
            raise ValueError("provide an optimization method")
        self.optimize_algo = optimize_algo
        
        # call functions
        self.variational_distribution
        self.sampling_latent
        self.generative_model
        self.optimize
        self.optimize_test
        self.imputation

    @define_scope
    def variational_distribution(self):
        """
        defines the variational distribution or inference network of the model
        q(z, l | x, s)


        """

        #q(z | x, s)
        if self.log_variational:
            x = tf.log(1 + self.expression)
        else:
            x = self.expression

        h = dense(x, self.n_hidden, activation=tf.nn.relu, \
                    bn=True, keep_prob=self.dropout_rate, phase=self.training_phase)
        for layer in range(2, self.n_layers + 1):
            h = dense(h, self.n_hidden, activation=tf.nn.relu, \
                bn=True, keep_prob=self.dropout_rate, phase=self.training_phase)

        
        self.qz_m = dense(h, self.n_latent, activation=None, \
                bn=False, keep_prob=None, phase=self.training_phase)
        self.qz_v = dense(h, self.n_latent, activation=tf.exp, \
                bn=False, keep_prob=None, phase=self.training_phase)
        
        # q(l | x, s)
        h = dense(x, self.n_hidden, activation=tf.nn.relu, \
                bn=True, keep_prob=self.dropout_rate, phase=self.training_phase)
        self.ql_m = dense(h, 1, activation=None, \
                bn=False, keep_prob=None, phase=self.training_phase)
        self.ql_v = dense(h, 1, activation=tf.exp, \
                bn=False, keep_prob=None, phase=self.training_phase)
    
    @define_scope
    def sampling_latent(self):
        """
        defines the sampling process on the latent space given the var distribution
        """
            
        self.z = gaussian_sample(self.qz_m, self.qz_v)
        self.library = gaussian_sample(self.ql_m, self.ql_v)
    
    @define_scope
    def generative_model(self):
        """
        defines the generative process given a latent variable (the conditional distribution)
        """
            
        # p(x | z, s)
        if self.batch is not None:
            h = tf.concat([self.z, self.batch], 1)
        else:
            h = self.z
        
        #h = dense(h, self.n_hidden,
        #          activation=tf.nn.relu, bn=True, keep_prob=self.dropout_rate, phase=self.training_phase)
        h = dense(h, self.n_hidden,
                  activation=tf.nn.relu, bn=True, keep_prob=None, phase=self.training_phase)
                
        for layer in range(2, self.n_layers + 1):
            if self.batch is not None:
                h = tf.concat([h, self.batch], 1)
            h = dense(h, self.n_hidden, activation=tf.nn.relu, \
                bn=True, keep_prob=self.dropout_rate, phase=self.training_phase)

        if self.batch is not None:
            h = tf.concat([h, self.batch], 1)        
        
        #mean gamma
        self.px_scale = dense(h, self.n_input, activation=tf.nn.softmax, \
                    bn=False, keep_prob=None, phase=self.training_phase)
        
        #dispersion
        if self.dispersion == "gene-cell":
            self.px_r = dense(h, self.n_input, activation=None, \
                    bn=False, keep_prob=None, phase=self.training_phase)
        elif self.dispersion == "gene":
            self.px_r = tf.Variable(tf.random_normal([self.n_input]), name="r")
        else:
            if self.batch_ind is None:
                raise ValueError("batch dispersion with no batch info")
            else:
                self.px_r = tf.Variable(tf.random_normal([self.num_batches, self.n_input]), name="r")

            
        #mean poisson
        self.px_rate = tf.exp(self.library) * self.px_scale

        #dropout
        if self.zi:
            self.px_dropout = dense(h, self.n_input, activation=None, \
                    bn=False, keep_prob=None, phase=self.training_phase)
        

    @define_scope
    def optimize(self):
        """
        write down the loss and the optimizer
        """
        
        # converting from batch to local quantities
        if self.dispersion == "gene-batch":
            local_dispersion = tf.matmul(self.batch, tf.exp(self.px_r))
        else: 
            local_dispersion = tf.exp(self.px_r)
            
        if self.library_mode == "numeric":
            local_l_mean = self.library_size_mean
            local_l_var = self.library_size_var
        else:
            local_l_mean = tf.matmul(self.batch, self.library_size_mean)
            local_l_var = tf.matmul(self.batch, self.library_size_var)
        
        
        # VAE loss
        if self.zi:
            recon = log_zinb_positive(self.expression, self.px_rate, local_dispersion, \
                                  self.px_dropout)
        else:
            recon = log_nb_positive(self.expression, self.px_rate, local_dispersion)
        
        kl_gauss_z = 0.5 * tf.reduce_sum(\
                        tf.square(self.qz_m) + self.qz_v - tf.log(1e-8 + self.qz_v) - 1, 1)
        kl_gauss_l = 0.5 * tf.reduce_sum(\
                        tf.square(self.ql_m - local_l_mean) / local_l_var  \
                            + self.ql_v / local_l_var \
                            + tf.log(1e-8 + local_l_var)  - tf.log(1e-8 + self.ql_v) - 1, 1)
        
        self.ELBO_gau = tf.reduce_mean(recon - self.kl_scale * kl_gauss_z - kl_gauss_l)
        
        # MMD loss
        if self.apply_mmd:
            self.mmd = mmd_objective(self.z, self.batch_ind, self.num_batches)
            self.loss = - self.ELBO_gau + self.mmd_scale *  self.mmd
        
        else:
            self.loss = - self.ELBO_gau
        
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        optimizer = self.optimize_algo
        with tf.control_dependencies(update_ops):
            self.train_step = optimizer.minimize(self.loss)
    
    @define_scope
    def optimize_test(self):
        # Test time optimizer to compare log-likelihood score of ZINB-WaVE
        update_ops_test = tf.get_collection(tf.GraphKeys.UPDATE_OPS, "variational")
        test_vars = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "variational")
        optimizer_test = tf.train.AdamOptimizer(learning_rate=0.001, epsilon=0.1)
        with tf.control_dependencies(update_ops_test):
            self.test_step = optimizer_test.minimize(self.loss, var_list=test_vars)
    
    @define_scope
    def imputation(self):
        # more information of zero probabilities
        if self.zi:
            self.zero_prob = tf.nn.softplus(- self.px_dropout + tf.exp(self.px_r) * self.px_r - tf.exp(self.px_r) \
                             * tf.log(tf.exp(self.px_r) + self.px_rate + 1e-8)) \
                             - tf.nn.softplus( - self.px_dropout)
            self.dropout_prob = - tf.nn.softplus( - self.px_dropout)

