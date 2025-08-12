
import sys
import os
import numpy as np
import pickle as pkl


file_path = sys.argv[1]
origin_tissue = sys.argv[2]
consensus_threshold = float(sys.argv[3])
outprefix = sys.argv[4]

with open(file_path, "rb") as file:
    met_times = pkl.load(file)

# get a dict mapping each migration to the time it occurred across all posterior samples
expected_migration_times_all = {}

primary_to_met = 0
met_to_met = 0
met_to_primary = 0
for graph in met_times.values():

    # only record binary classifications per graph sample
    pm = False
    mm = False
    mp = False

    for migration, time in graph.items():

        middle_time = (time[0] + time[1]) / 2

        # track expected migration times
        if migration not in expected_migration_times_all:
            expected_migration_times_all[migration] = []
        expected_migration_times_all[migration].append(middle_time)

        # get classifications
        source, recipient = migration.split("_")[0:2]
        if source == origin_tissue and pm == False:
            primary_to_met += 1
            pm = True
        elif recipient == origin_tissue and mp == False:
            met_to_primary += 1
            mp = True
        elif source != origin_tissue and recipient != origin_tissue and mm == False:
            met_to_met += 1
            mm = True

# find only the edges above the consensus threshold and get the expected time for each migration
total_graphs = len(met_times)
expected_migration_times_all = {
    k: v
    for k, v in expected_migration_times_all.items()
    if len(v) / total_graphs > consensus_threshold
}

# Merge edges for the same source and recipient to ignore the multiedge aspect
expected_migration_times_no_multiedge = {}
for migration, time in expected_migration_times_all.items():
    no_multiedge = "_".join(migration.split("_")[0:2])
    if no_multiedge not in expected_migration_times_no_multiedge:
        expected_migration_times_no_multiedge[no_multiedge] = []
    expected_migration_times_no_multiedge[no_multiedge].extend(time)


expected_migration_times_all = {
    k: np.mean(v) for k, v in expected_migration_times_all.items()
}
expected_migration_times_no_multiedge = {
    k: np.mean(v) for k, v in expected_migration_times_no_multiedge.items()
}

# get the consensus classifications
p_to_m = (primary_to_met / total_graphs) > consensus_threshold
m_to_m = (met_to_met / total_graphs) > consensus_threshold
m_to_p = (met_to_primary / total_graphs) > consensus_threshold

# write out the expected migration times
consensus_threshold = int(consensus_threshold * 100)
outname = os.path.basename(outprefix)
with open(
    outprefix + f"_expected_migration_times_no_multiedge_{consensus_threshold}.csv", "w"
) as file:
    file.write("name,source_recipient,mid_time\n")
    for migration, time in expected_migration_times_no_multiedge.items():
        file.write(f"{outname},{migration},{time}\n")

with open(
    outprefix + f"_expected_migration_times_all_{consensus_threshold}.csv", "w"
) as file:
    file.write("name,source_recipient_multiedgeNum,mid_time\n")
    for migration, time in expected_migration_times_all.items():
        file.write(f"{outname},{migration},{time}\n")

# write out the consensus classifications
with open(
    outprefix + f"_consensus_classifications_{consensus_threshold}.csv", "w"
) as file:
    file.write("name,met_to_met,met_to_primary\n")
    file.write(f"{outname},{m_to_m},{m_to_p}")
