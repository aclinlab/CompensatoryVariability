#!/usr/bin/env python

import sys
import os

import h5py
import neuprint
from neuprint import Client, queries, SegmentCriteria
from extract_labels_from_volume import extract_labels_from_volume

import numpy as np

c = Client('neuprint.janelia.org', dataset='hemibrain:v1.1', token='eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImFuZHJldy5saW5Ac2hlZmZpZWxkLmFjLnVrIiwibGV2ZWwiOiJub2F1dGgiLCJpbWFnZS11cmwiOiJodHRwczovL2xoMy5nb29nbGV1c2VyY29udGVudC5jb20vYS0vQUF1RTdtQmp2LXlGaUlDWjdFQ3AwUHhzYmNVSHNDclppTU96eWtHbnFSemk_c3o9NTA_c3o9NTAiLCJleHAiOjE3NjE0MzE3Mzh9.kfurqTReeqmGsFHQdkKKGhTcaw3JR-FnpxCpMpS7uOA')
print(c.fetch_version())

print("Loading volume")
with h5py.File('hemibrain-v1.1-primary-roi-segmentation.h5', 'r') as f:
    roi_vol_scale_5 = f['volume-256nm'][:]

KCs = SegmentCriteria(type="^KC.*",regex=True)
KCneurons, KCroicounts = neuprint.queries.fetch_neurons(KCs)
KCids = KCneurons['bodyId']
KCidsarray = np.asarray(KCids)
np.savetxt("KCids_v1.1.csv",KCidsarray,delimiter=",",fmt="%u")
for i in range(0,len(KCids)):
    print(f"Fetching skeleton for KC {i}")
    KC = SegmentCriteria(bodyId=KCids[i])
    KCskel=c.fetch_skeleton(KCids[i],format='pandas',heal=True)
    print(f"Extracting ROI for each skeleton point for KC {i}")
    extract_labels_from_volume(KCskel, roi_vol_scale_5, vol_scale=5, label_names=c.primary_rois)
    KCskel.rename(columns={'label_name': 'roi'}, inplace=True)
    KCskel.to_csv(path_or_buf="KCskel"+str(i)+".csv")
    KCinputs = neuprint.queries.fetch_synapse_connections(target_criteria=KC)
    print(f"Found {len(KCinputs)} inputs to KC {i} bodyId {KCids[i]}")
    KCinputs.to_csv(path_or_buf="KCinputs"+str(i)+".csv")
    KCoutputs = neuprint.queries.fetch_synapse_connections(source_criteria=KC)
    print(f"Found {len(KCoutputs)} outputs from KC {i} bodyId {KCids[i]}")
    KCoutputs.to_csv(path_or_buf="KCoutputs"+str(i)+".csv")
