# Copyright 2015 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Example / benchmark for building a PTB LSTM model.

The hyperparameters used in the model:
- init_scale - the initial scale of the weights
- learning_rate - the initial value of the learning rate
- max_grad_norm - the maximum permissible norm of the gradient
- num_layers - the number of LSTM layers
- num_steps - the number of unrolled steps of LSTM
- hidden_size - the number of LSTM units
- max_epoch - the number of epochs trained with the initial learning rate
- max_max_epoch - the total number of epochs for training
- keep_prob - the probability of keeping weights in the dropout layer
- lr_decay - the decay of the learning rate for each epoch after "max_epoch"
- batch_size - the batch size

To run:

$ python3 lstm.py --data_path=simple-examples/data/ --config_file=config.txt --checkpoint_path=path/

To install tensorflow from source:
run in the root of the tree:

$ ./configure
/Library/Frameworks/Python.framework/Versions/3.4/bin/python3

$ bazel build -c opt //tensorflow/tools/pip_package:build_pip_package
# ...or, with GPU support
$ bazel build -c opt --config=cuda //tensorflow/tools/pip_package:build_pip_package

$ bazel-bin/tensorflow/tools/pip_package/build_pip_package /tmp/tensorflow_pkg

# The name of the .whl file will depend on your platform.
$ sudo pip3 install /tmp/tensorflow_pkg/tensorflow-0.9.0-py3-none-any.whl

Seg fault during tensorflow import:
Having the same problem, It appears that it fails when loading libcuda.dylib. Unfortunately, since #2878 it's now looking for libcuda.1.dylib instead.
Creating a symlink or changing the name of libcuda.dylib to libcuda.1.dylib in /usr/local/cuda/lib seems to do the job.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time
from random import sample
import os, sys, glob

import numpy as np
import tensorflow as tf

from tensorflow.python.ops import rnn_cell
from tensorflow.python.ops import rnn

flags = tf.flags
logging = tf.logging

flags.DEFINE_string(
		"model", "small",
		"A type of model. Possible options are: small, medium, large.")
flags.DEFINE_string("data_path", None, "data_path")
flags.DEFINE_string("config_file", None, "config_file")
flags.DEFINE_string("checkpoint_path", None, "checkpoint_path")
flags.DEFINE_string("result_path", None, "result_path")

FLAGS = flags.FLAGS

class LSTM(object):

	def __init__(self, is_training, config):
		self.batch_size = batch_size = config.batch_size
		self.num_steps = num_steps = config.num_steps
		self.leave_out = leave_out = config.leave_out
		size = config.hidden_size
		self.sample_dim = sample_dim = config.sample_dim

		self._input_data = inputs = tf.placeholder(tf.float32, [batch_size, num_steps, sample_dim])
		self._targets = tf.placeholder(tf.float32, [batch_size, num_steps, sample_dim])
		self._seq_l_matrix = tf.placeholder(tf.bool, [batch_size, num_steps, sample_dim])
		self._seq_length = tf.placeholder(tf.int32, [batch_size])

		# Slightly better results can be obtained with forget gate biases
		# initialized to 1 but the hyperparameters of the model would need to be
		# different than reported in the paper.
		lstm_cell = tf.nn.rnn_cell.BasicLSTMCell(size, forget_bias=config.forget_bias)
		if is_training and config.keep_prob < 1:
			lstm_cell = tf.nn.rnn_cell.DropoutWrapper(
					lstm_cell, output_keep_prob=config.keep_prob)
		cell = rnn_cell.MultiRNNCell([lstm_cell] * config.num_layers)

		self._initial_state = cell.zero_state(batch_size, tf.float32)

		inputs = [tf.squeeze(input_, [1])
			for input_ in tf.split(1, num_steps, inputs)]

		targets = self._targets

		# print("Inputs")
		# print(inputs)
		# print("Targets")
		# print(targets)


		softmax_w = tf.get_variable("RNN/softmax_w", [sample_dim, size])
		softmax_b = tf.get_variable("RNN/softmax_b", [1,sample_dim])

		def transform_output_to_input(input_, output, time):
			output_t = tf.transpose(output)

			softmax_w = tf.get_variable("softmax_w", [sample_dim, size])
			softmax_b = tf.get_variable("softmax_b", [1,sample_dim])

			new_input = tf.transpose(tf.matmul(softmax_w, output_t)) + softmax_b

			return(new_input)

		outputs, state = rnn.rnn(cell, inputs, initial_state=self._initial_state, 
						sequence_length=self._seq_length, feed_previous = True,
						function_on_input = transform_output_to_input)

		# print("outputs")
		# print(outputs[0])

		outputs = [tf.transpose(o) for o in outputs]

		predictions = tf.pack(outputs)

		# print("Predictions")
		# print(predictions)

		predictions = tf.reshape(tf.transpose(predictions), [batch_size, size, num_steps])

		predictions = tf.unpack(predictions)

		# for each sample in batch
		predictions = [tf.transpose(tf.matmul(softmax_w, one_prediction)) + tf.gather(softmax_b, [0] * num_steps) 
							for one_prediction in predictions]

		# print("bunny")
		# print(predictions)

		predictions = tf.pack(predictions)

		loss = tf.square(predictions - targets)
		loss = tf.select(self._seq_l_matrix, loss, tf.zeros([batch_size, num_steps, sample_dim], tf.float32))

		self._cost = cost = tf.log(tf.reduce_sum(loss) / batch_size / sample_dim)
		self._final_state = state
		self._loss = tf.reduce_sum(loss)
		self._predictions = predictions
		self.global_step = tf.Variable(0, name='global_step', trainable=False)

		if not is_training:
			return

		self._lr = tf.Variable(0.0, trainable=False)
		#tvars = tf.trainable_variables()
		#grads = tf.gradients(cost, tvars)
		
		self._train_op = tf.train.AdamOptimizer(self.lr).minimize(cost)
		#self._train_op = optimizer.apply_gradients(zip(grads, tvars))

	def assign_lr(self, session, lr_value):
		session.run(tf.assign(self.lr, lr_value))

	def assign_global_step(self, session, step):
		session.run(tf.assign(self.global_step, step))

	@property
	def input_data(self):
		return self._input_data

	@property
	def targets(self):
		return self._targets

	@property
	def initial_state(self):
		return self._initial_state

	@property
	def seq_l_matrix(self):
		return self._seq_l_matrix

	@property
	def seq_length(self):
		return self._seq_length

	@property
	def loss(self):
		return self._loss

	@property
	def predictions(self):
		return self._predictions

	@property
	def loss(self):
		return self._loss

	@property
	def cost(self):
		return self._cost

	@property
	def final_state(self):
		return self._final_state

	@property
	def lr(self):
		return self._lr

	@property
	def train_op(self):
		return self._train_op


class config_reader:
	def __init__(self, config_file):
		with open(config_file) as f:
			for l in f:
				tokens = l.strip().split(" = ")
				attr_name, value = tokens
				if (value.find('.') == -1):
					value = int(value)
				else:
					value = float(value)
				setattr(self, attr_name, value)

class reader(object):
	@staticmethod
	def split_into_batches(data, data_lengths, batch_size, num_steps, leave_out, sample_dim):
		l_data = len(data)
		per_batch=[]

		for i in range(l_data // batch_size):
			x = data[(i*batch_size) : ((i+1)*batch_size)]
			lengths = data_lengths[(i*batch_size) : ((i+1)*batch_size)]

			# Do not look at predictions at first leave_out time points and after the last time point of the series
			l_matrix = []
			for l in lengths:
				new = []
				new.extend([[False] * sample_dim for j in range(leave_out)])
				new.extend([[True] * sample_dim for j in range(leave_out, l-1)])
				new.extend([[False] * sample_dim for k in range(max(l-1,leave_out), num_steps)])
				l_matrix.append(new)

			# switching the last two dimentions of x
			new_x = []
			for sample in x:
				new_sample = []
				for j in range(num_steps):
					new_sample.append([trajectory[j] for trajectory in sample])

				new_x.append(new_sample)
			x = new_x


			# Predict the next element
			y = []
			for sample in x:
				new_label = sample[1:]
				new_label.append([0] * sample_dim)
				y.append(new_label)

			x = [np.reshape(sample, [num_steps, sample_dim]) for sample in x]
			y = [np.reshape(sample, [num_steps, sample_dim]) for sample in y]
			l_matrix = [np.reshape(sample, [num_steps, sample_dim]) for sample in l_matrix]

			per_batch.append((x, y, lengths, l_matrix))

		return per_batch

	@staticmethod
	def raw_data(data_path, percent_train, percent_valid, sample_dim):
		data=[]
		sample_list = []

		for data_file in os.listdir(data_path):
			data_file = os.path.join(data_path, data_file)
			if not data_file.endswith(".txt"):
				next

			data_item = []
			n_trajectories = 0
			with open(data_file) as f:
				for l in f:
					data_item.append(list(map(float,l.strip().split("\t"))))
					n_trajectories += 1
					if (n_trajectories == sample_dim):
						break
			sample_list.append(data_file)

			#assert(len(data_item) == sample_dim)
			data.append(data_item)

		# number of samples
		l_data = len(data)
		lengths = [len(sample[0]) for sample in data]

		train_indices = list(range(int(l_data * percent_train)))
		valid_indices = list(range(int(l_data * percent_train), int(l_data * (percent_train + percent_valid))))
		test_indices = list(range(int(l_data * (percent_train + percent_valid)), l_data))

		train_data = [data[i] for i in train_indices]
		valid_data = [data[i] for i in valid_indices]
		test_data = [data[i] for i in test_indices]

		train_l = [lengths[i] for i in train_indices]
		valid_l = [lengths[i] for i in valid_indices]
		test_l = [lengths[i] for i in test_indices]

		train_sample_list = [sample_list[i] for i in train_indices]
		valid_sample_list = [sample_list[i] for i in valid_indices]
		test_sample_list = [sample_list[i] for i in test_indices]

		return [train_data, valid_data, test_data], [train_l, valid_l, test_l], [train_sample_list, valid_sample_list, test_sample_list]

	@staticmethod
	def padding(data, num_steps):
		for i in range(len(data)):
			for j in range(len(data[i])):
				data[i][j].extend([0] * (num_steps - len(data[i][j])))

		return(data)

def save_test_predictions(test_data, test_sample_names, predictions, path = "predictions_test_set"):
	if not os.path.exists(path):
		os.makedirs(path)

	data_path = os.path.join(path, "initial_data", "")
	if not os.path.exists(data_path):
		os.makedirs(data_path)
	prediction_path = os.path.join(path, "predictions", "")
	if not os.path.exists(prediction_path):
		os.makedirs(prediction_path)

	for i, sample_name in enumerate(test_sample_names):
		with open(data_path + os.path.basename(sample_name), "w") as f:
			for trajectory in test_data[i]:
				f.write("\t".join([str(trajectory[j]) for j in range(len(trajectory))]) + "\n")

		with open(prediction_path + os.path.basename(sample_name), "w") as f:
			for trajectory in np.transpose(predictions[i]):
				f.write("\t".join([str(trajectory[j]) for j in range(len(trajectory))]) + "\n")


def run_epoch(session, m, data, data_lengths, eval_op, verbose=False):
	"""Runs the model on the given data."""
	epoch_size = len(data) // m.batch_size
	start_time = time.time()
	costs = 0.0
	iters = 0
	predictions = []

	state = m.initial_state.eval()
	for step, (x, y, l, l_matrix) in enumerate(
		reader.split_into_batches(data, data_lengths, m.batch_size, m.num_steps, m.leave_out, m.sample_dim)):
		
		cost, final_state, _, loss, prediction = session.run([m.cost, m.final_state, eval_op, m.loss, m.predictions],
																 {m.input_data: x,
																	m.targets: y,
																	m.initial_state: state,
																	m.seq_l_matrix: l_matrix, 
																	m.seq_length: l})
		costs += cost
		iters += 1
		predictions.append(prediction[0])

		if verbose and (iters < 10 or (loss == 0)):
			print(iters)
			print("%.3f loss: %.3f speed: %.0f wps" %
						(step * 1.0 / epoch_size, costs / iters,
						 iters * m.batch_size / (time.time() - start_time)))
			print("Cost")
			print(cost)
			print(loss)
			print("predictions")
			print(predictions)
			print("labels")
			print(y)
			print("initial_data")
			print(x)
			print(l)
			print(sum(l))

	return costs / iters, predictions


def main(_):
	if not FLAGS.data_path:
		raise ValueError("Must set --data_path data directory")

	if not FLAGS.config_file:
		raise ValueError("Must set --config_file config file")

	checkpoint_path = FLAGS.checkpoint_path

	if (checkpoint_path):
		config_file_name = os.path.basename(FLAGS.config_file)
		if config_file_name.endswith('.txt'):
			config_file_name = config_file_name[:-4]

		script_name = os.path.basename(__file__)[:-3]

		checkpoint_path = os.path.join(checkpoint_path, config_file_name, "")
		if not os.path.exists(checkpoint_path):
			os.makedirs(checkpoint_path)

	result_path = "predictions_test_set"
	if FLAGS.result_path:
		result_path = FLAGS.result_path
		
	config = config_reader(FLAGS.config_file)
	eval_config = config_reader(FLAGS.config_file)
	eval_config.batch_size = 1
	eval_config.max_max_epoch = 1

	raw_data, data_lengths, samples = reader.raw_data(FLAGS.data_path, 0.6, 0.2, config.sample_dim)
	train_data, valid_data, test_data = raw_data
	train_lengths, valid_lengths, test_lengths = data_lengths
	train_samples, valid_samples, test_samples = samples

	with tf.Graph().as_default(), tf.Session() as session:
		initializer = tf.random_uniform_initializer(-config.init_scale, config.init_scale)
		
		with tf.variable_scope("model", reuse=None, initializer=initializer):
			m = LSTM(is_training=True, config=config)
		with tf.variable_scope("model", reuse=True, initializer=initializer):
			mvalid = LSTM(is_training=False, config=config)
			mtest = LSTM(is_training=False, config=eval_config)

		train_data = reader.padding(train_data, m.num_steps)
		valid_data = reader.padding(valid_data, mvalid.num_steps)
		test_data = reader.padding(test_data, mtest.num_steps)

		tf.initialize_all_variables().run()
		saver = tf.train.Saver()

		if (checkpoint_path):
			ckpt = tf.train.get_checkpoint_state(checkpoint_path)
			if ckpt and ckpt.model_checkpoint_path:
				saver.restore(session, ckpt.model_checkpoint_path)
				print("Model restored from %s Step step: %d" % (checkpoint_path, m.global_step.eval()))

		lr = config.learning_rate
		for i in range(m.global_step.eval(), config.max_max_epoch):
			lr_decay = config.lr_decay ** max(i - config.max_epoch, 0.0)
			if ( lr  > config.lr_stop_decay):
				lr =  config.learning_rate * lr_decay
			m.assign_lr(session, lr)

			#print("Epoch: %d Learning rate: %.3f" % (i + 1, session.run(m.lr)))
			assert(len(train_data) == len(train_lengths))

			shuffle_indices  = list(range(len(train_data)))
			shuffle_indices = sample(shuffle_indices, len(shuffle_indices))
			train_data_shuffled = [train_data[j] for j in shuffle_indices]
			train_lengths_shuffled = [train_lengths[j] for j in shuffle_indices]

			train_loss, _ = run_epoch(session, m, train_data_shuffled, train_lengths_shuffled, m.train_op, verbose=False)
			#print("Epoch: %d Train loss: %.3f" % (i + 1, train_loss))
			valid_loss, _ = run_epoch(session, mvalid, valid_data, valid_lengths, tf.no_op())
			print("Epoch: %d Learning rate: %.3f Train loss: %.3f Valid loss: %.3f" % 
				(i + 1, session.run(m.lr), train_loss, valid_loss))

			m.assign_global_step(session, i+1)

			if (checkpoint_path and (i+1) % 5 == 0):
				save_path = saver.save(session, checkpoint_path, global_step=i+1)
				print("Model saved in file: %s. Step: %d" % (save_path, m.global_step.eval()))


		test_loss, predictions = run_epoch(session, mtest, test_data, test_lengths, tf.no_op(), verbose=False)
		save_test_predictions(test_data, test_samples, predictions, result_path)
		print("Test loss: %.3f" % test_loss)


if __name__ == "__main__":
	tf.app.run()
