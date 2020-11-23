from wgs.workflows.cohort_qc import tasks
import sys

mafs = sys.argv[1]
mafs = list(map(str, mafs.strip('[]').split(',')))
mafs = [m.strip().replace("'", "") for m in mafs]

labels = sys.argv[2]
#labels = [labels] * len(mafs)
labels = list(map(str, labels.strip('[]').split(',')))
labels = [m.strip().replace("'", "") for m in labels]


out = sys.argv[3]

tasks.merge_mafs(mafs, out, labels)
