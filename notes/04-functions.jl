using Resonance

# - Takes ~ 6min for 1500 samples:
# - Writes to ENV["SCRATCH_SPACE"], which is "./scratch/" unless otherwise set
Resonance.write_gfs_arrow()
Resonance.write_gfs_arrow(; kind="ecs_rename", stratified=true)
Resonance.write_gfs_arrow(; kind="kos_rename", stratified=true)
Resonance.write_gfs_arrow(; kind="pfams_rename", stratified=true)
