#!/usr/bin/env python3

import baltic as bt
import sys
import datetime

tree_path = sys.argv[1]

tree = bt.loadNexus(tree_path)
tree.sortBranches()

leaf_node_times = [i.absoluteTime for i in tree.Objects if i.branchType == "leaf" and i.traits["location"] != "Outside_San_Diego"]

mrca = tree.commonAncestor([i for i in tree.Objects if i.branchType == "leaf" and i.traits["location"] != "Outside_San_Diego"])

time_interval = [
    mrca.absoluteTime,
    max(leaf_node_times)
]

def get_iso_week_start(x):
    calendar_date = [int(i) for i in bt.calendarDate(x).split("-")]
    iso_week = datetime.date(*calendar_date).isocalendar()[1]
    start_date_iso_week = datetime.datetime.fromisocalendar(calendar_date[0], iso_week, 1).strftime("%Y-%m-%d")
    return bt.decimalDate(start_date_iso_week)

iso_weeks = [get_iso_week_start(time_interval[0]), get_iso_week_start(time_interval[1])]

evaluation_times = []
end = iso_weeks[1]
while end > iso_weeks[0]:
    evaluation_times.append(end)
    end -= (1/52)

evaluation_times = evaluation_times[::-1]

ancestral_times = [time_interval[1]-i for i in evaluation_times[:-1]]
evaluation_times = [time_interval[1]-i for i in  evaluation_times[1:]]

print_times = lambda x: ",".join([str(i) for i in x])

print("-evaluationTime {} -ancestralTime {}".format(print_times(evaluation_times), print_times(ancestral_times)))
