import re
import math
import time
import os
import glob
import os.path
import logging
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import scipy.stats
import h5py
import matplotlib.pyplot as plt
import quickpen
import penplot
import time
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
    info["TrajectoryLines"]=include_trajectory_lines()
    info["EndTime"]=include_end_time()
    info["TotalInfected"]=include_total_infected()

    text="""\\documentclass{{article}}
\\usepackage{{graphicx}}
\\usepackage{{hyperref}}
\\begin{{document}}
\\title{{Report on {Title}}}
\\date{{\\today}}
\\author{{Generated with \\texttt{{report.py}}}}
\\maketitle

{CommandlineOptionsTable}

{CodeTraitsTable}

{TrajectoryLines}

{EndTime}

{TotalInfected}

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

def include_end_time():
    return include_figure("end_time_hist.pdf", "0.7",
        "Each bar shows a count of how many realizations completed "+
        "at a given time.", "fig:endtimes")

def include_total_infected():
    return include_figure("total_infected_hist.pdf", "0.7",
        "Each bar shows the total number of realizations whose "
        "infections fell in the given range", "fig:totalinfected")

def trajectory_lines(h5f):
    penplot.plot_trajectory_lines(quickpen.FileTrajectories(h5f))

def include_trajectory_lines():
    return include_figure("trajectory_lines.pdf", "0.7",
        "Each line represents a separate realization from the simulation."+
        "This shows all infected individuals, whether exposed or infectious.",
        "fig:trajlines")

def include_figure(fileglob, scale, caption, label):
    snippet=""
    targets=glob.glob(fileglob)
    blanks={"scale" : scale, "caption" : caption}
    for idx, f in enumerate(targets):
        blanks["filename"]=f
        if len(targets)>1:
            blanks["label"]="{0}{1}".format(label, idx)
        else:
            blanks["label"]=label
        snippet=snippet+"""\\begin{{figure}}
\\centerline{{\\includegraphics[scale={scale}]{{{filename}}}}}
\\caption{{{caption}\\label{{{label}}}}}
\\end{{figure}}
""".format(**blanks)

    return snippet


def make_report(filename):
    f=h5py.File(filename, "r")
    info=dict()
    info.update(metadata(f))
    info["Title"]="{0}--{1}".format(filename, info["UUID"][0:5])
    outfile="report.tex"
    write_report(info, outfile)


def parallel_generate(filename, timeout=10):
    '''
    Timeout in minutes.
    '''
    todo=[
      ["python", "report.py", "--file", filename, "--multiples"],
      ["python", "report.py", "--file", filename, "--summary"],
      ["python", "report.py", "--file", filename, "--lines"]
    ]
    processes=list()
    for t in todo:
        processes.append(subprocess.Popen(t))
    interval=2 # seconds
    for check in range(int(timeout*60/interval)):
        done=list()
        for idx, proc in enumerate(processes):
            if proc.poll()==0:
                done.append(idx)
                logging.info("{0} completed".format(todo[idx][4]))
        done.reverse()
        for d in done:
            processes.pop(d)
            todo.pop(d)
        if len(processes)==0:
            return True
        time.sleep(interval)
    return False


if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO)
    parser=DefaultArgumentParser(description="Quick look at an H5 file")
    parser.add_argument("--file", dest="file", action="store",
        default="rider.h5", help="data file to read")
    parser.add_function("report", "Build the report")
    parser.add_function("lines", "Write the image of all realizations")
    parser.add_function("multiples", "Small multiples graph of one realization")
    parser.add_function("summary", "Several summary graphs")
    parser.add_function("generate", "Create graphs and report")

    args=parser.parse_args()

    filename=args.file
    if args.report:
        make_report(filename)
        # if not parser.any_function():
        #     parser.print_help()
    elif args.lines:
        f=h5py.File(filename, "r")
        trajectory_lines(f)
    elif args.multiples:
        f=h5py.File(filename, "r")
        single_trajectory_small_multiples(f)
    elif args.summary:
        f=h5py.File(filename, "r")
        summaries(f)
    elif args.generate:
        parallel_generate(filename)

