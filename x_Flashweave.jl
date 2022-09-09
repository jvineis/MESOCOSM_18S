#!/usr/bin/env julia

using FlashWeave
println("Hello climber, I'm running Flashweave for you")

v4_count_data = "x_SWARM-counts-for-flashweave-with-MAGs-transposed-min50.csv"
v4_meta = "x_SWARM-counts-for-flashweave-transposed-min50-metadata.csv"
v4_network_count = learn_network(v4_count_data,v4_meta,n_obs_min=15)
save_network("x_SWARM-counts-for-flashweave-with-MAGs-output.edgelist", v4_network_count)
save_network("x_SWARM-counts-for-flashweave-with-MAGs-output.gml", v4_network_count)