#!/usr/bin/env bash

# run mongo daemon in the background
mongod --dbpath . --port 1234 &

# FIXME: find a way to get the number of available gpus and their ids
# launch a hyperopt worker on each gpu
for  gpu in gpus
do
	CUDA_VISIBLE_DEVICE=gpu PYTHONPATH=.. hyperopt-mongo-worker --mongo=localhost:1234/scvi_db --poll-interval=0.1 &
done

