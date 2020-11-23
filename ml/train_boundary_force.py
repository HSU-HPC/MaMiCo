import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import numpy as np
from glob import glob

#dataset = tf.data.experimental.make_csv_dataset('argon.boundary.samples.txt', batch_size=32, field_delim=';',
#            column_names=['distance', 'outer', 'inner', 'total'], select_columns=['distance', 'outer'], num_epochs=10)

filenames = sorted(glob('md30_data/*csv'))
filenames = tf.convert_to_tensor(filenames, dtype=tf.string)

dataset = tf.data.experimental.CsvDataset(filenames, 
            [tf.float64 for x in range(18)], field_delim=',').shuffle(1024)

def preprocess(*fields):
    return tf.stack([fields[0], fields[10], fields[11], fields[12]]), tf.stack(fields[-3:])

train_ds = dataset.map(preprocess).take(100000).batch(128).repeat(1)
val_ds = dataset.map(preprocess).skip(1000).take(30000).batch(128).repeat(1)
test_ds = dataset.map(preprocess).skip(1000).take(10000).batch(1000).repeat(1)


inputs = layers.Input(shape=(4,), name='input')
normalize = layers.BatchNormalization(center=False, scale=False)(inputs)
dense = layers.Dense(24, activation="tanh")(normalize)
dense = layers.Dense(10, activation="tanh")(dense)
outputs = layers.Dense(3, activation=None, name='output')(dense)


model = keras.Model(inputs=inputs, outputs=outputs, name="boundary_force_model")

model.summary()

model.compile('Adam', keras.losses.mean_squared_error)

model.fit(train_ds, validation_data=val_ds, epochs=1, callbacks=[
    keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
])


result = model.evaluate(test_ds)

model.save("./saved_model/", include_optimizer=False, save_format='tf')