# SangerBlobML
Tensorflow based model to identify dye blobs in sanger trace data. 

Usage: python SnagerBlobML.py

Model: "sequential_44"
_________________________________________________________________
 Layer (type)                Output Shape              Param #
=================================================================
 layer_normalization_25 (Lay  (None, 5, 40)            80
 erNormalization)

 flatten_30 (Flatten)        (None, 200)               0

 dense_61 (Dense)            (None, 32)                6432

 dense_62 (Dense)            (None, 1)                 33

=================================================================
Total params: 6,545
Trainable params: 6,545
Non-trainable params: 0
_________________________________________________________________
