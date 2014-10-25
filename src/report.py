import re
import math
import logging
import xml.etree.ElementTree as ET
import numpy as np
import scipy.stats
import h5py
import matplotlib.pyplot as plt
import quickpen
import penplot
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)

def metadata(h5f):
    info={}
    attrs=h5f['/trajectory'].attrs
    info['CompileTime']=attrs['COMPILETIME'].tostring().decode("utf-8")
    info['Makefile']=attrs['CONFIG'].tostring().decode("utf-8")
    gitrepo=attrs['VERSION'].tostring().decode("utf-8")
    m=re.search("@([a-zA-Z0-9\-\.]+):(\w+)\/(\w+)\.git:(\w+)",gitrepo)
    if m:
        website, who, project, version=m.groups()
        info['RepoUrl']='https://{0}/{1}/{2}/commit/{3}'.format(website, who, project, version)
        info['RepoEnglish']="{0}/{1} project at {2}".format(who, project, website)
    s, e, i, r=attrs['Initial Values']
    info['InitialValues']='Susceptible {0} exposed {1} infectious {2}  recovered {3}'.format(
        s, e, i, r)
    root=ET.fromstring(attrs['Options'].tostring().decode("utf-8"))
    commandline_options=dict()
    for option in root.iter("option"):
        name=option.find("name").text
        values=list()
        for v in option.find("values").findall("value"):
            values.append(v.text)
        commandline_options[name]=", ".join(values)

    info['CommandlineOptions']=commandline_options
    print(type(attrs['uuid']))
    print(attrs['uuid'])
    print(attrs['uuid'].shape)
    info['UUID']=attrs['uuid'].tostring().decode("utf-8")
    return info


def write_report(info, outfile):
    options_table=["\\begin{tabular}{ll}", "Option & Value \\\\ \\hline"]
    for k, v in sorted(info["CommandlineOptions"].items()):
        options_table.append("{0} & {1} \\\\".format(k, v))
    options_table.append("\\end{tabular}")
    options_table="\n".join(options_table)
    info["CommandlineOptionsTable"]=options_table

    provenance_table=["\\begin{tabular}{ll}", "Trait & Value \\\\ \\hline"]
    provenance_table.append("Compile time & {0} \\\\".format(info["CompileTime"]))
    provenance_table.append("Code repository & \\href{{{0}}}{{{1}}} \\\\".format(info["RepoUrl"], info["RepoEnglish"]))
    provenance_table.append("Initial values & {0} \\\\".format(info["InitialValues"]))
    provenance_table.append("Unique Tag & {0} \\\\".format(info["UUID"]))
    provenance_table.append("\\end{tabular}")
    info["CodeTraitsTable"]="\n".join(provenance_table)

    text="""\\documentclass{{article}}
\\usepackage{{hyperref}}
\\begin{{document}}
\\title{{Report on {Title}}}
\\date{{\\today}}
\\author{{Generated with \\texttt{{report.py}}}}
\\maketitle

{CommandlineOptionsTable}

{CodeTraitsTable}

\\end{{document}}
""".format(**info)
    f=open(outfile, "w")
    f.write(text)


def single_trajectory_small_multiples(h5f):
    dsetname=quickpen.trajectories(h5f)[0]
    tr,times,obs=quickpen.per_pen_trajectory(h5f, dsetname)
    penplot.small_multiples(tr, times, obs, 4)

def summaries(h5f):
    aggregate=quickpen.summary_of_ensemble(h5f, 50)
    penplot.total_infected_plot(aggregate[0].data)
    penplot.end_time_plot(aggregate[1].data)
    all_states_obj=aggregate[2]
    print(type(all_states_obj))
    seir=all_states_obj.seir()
    logger.debug("seir shape {0}".format(seir.shape))
    times=all_states_obj.times()
    logger.debug("times shape {0}".format(times.shape))
    penplot.trajectory_density_plot(seir[:,1], times, "Exposed")
    penplot.trajectory_density_plot(seir[:,2], times, "Infected")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser=DefaultArgumentParser(description="Quick look at an H5 file")
    parser.add_argument("--file", dest="file", action="store",
        default="rider.h5", help="data file to read")

    args=parser.parse_args()

    filename=args.file
    f=h5py.File(filename, "r")
    info=dict()
    info.update(metadata(f))
    info["Title"]="{0}--{1}".format(filename, info["UUID"][0:5])
    outfile="report.tex"
    write_report(info, outfile)
    # if not parser.any_function():
    #     parser.print_help()
