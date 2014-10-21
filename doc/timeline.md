# Timeline for Development of Code

## Monday 13 October 2014

Today is the first day to write code.
I've spent plenty of time thinking about the model. There is
a description of a feedlot that's been looked-over by a few people.
The plan is to use a different project,
[SIRS Demonstration](https://github.com/adolgert/sirs), as a
boilerplate. I will make a series of models which can then be
used for comparison.

Normally the first steps would be to make a simple model and
validate it. In this case, I have some validation on the SIRS
model already, which is close to the basic model, and the big
question about modeling with this technique is whether it will
scale to a significant number of cattle. The model which is most
simple theoretically, an all-to-all well-mixed infection model,
is the slowest version. For N cattle, there are N^2 transitions
with this model. The plan is to push, as quickly as possible, to
a model that is more scalable, and then go back and look more
closely at the transition distributions.

### Copy SIRS Demo code

It turns out that the SIRS demo code needed some tweaking.
Once that runs, I have to turn it into an SEIR code. The
pathology of FMDV in an animal looks like SEIR for the basic
question of infectiousness. Whether the animal is clinical
is a second question, which we'll skip right now.
Plus, let's use exponential distributions for everything
and go back later to put in the more detailed distributions,
because I'm not sure the Gamma distribution is implemented
in the current code base.

There is always this moment, when changing a four hundred
line C++ code, where I wonder whether it will ever compile again.

### Timing with All-to-all

This SEIR model is a basic all-to-all, and it runs 1024 animals
in about 30s. That's about what I expected. The executable
is called "seir".


## Tuesday 14 October

### Adding Pens
The first step to scalability is to treat within-pen transitions
separately from pen-to-pen transitions. Let's modify the all-to-all
so that, during construction of the net, the transitions between
individuals are different depending on their pen-to-pen distance.

This runs, but no faster, of course.

### Coalescing Transitions
I expect some improvement in speed from a simple trick. Each
individual will have its set of compartments, SEIR, where
the E-I and I-R transitions are non-exponential, determined
by Mardones et al. We can, however, treat the infection process
in a grouped way. Each infectious individual is in contact with
many susceptibles, and we can draw one random number to
determine whether it infects any of them, then a second random
number to determine which is fired. This is different from
the normal individual-to-individual transition where choosing
to infect predetermines who gets infected.

Took some time to get this technique right. It seems slow, taking
more than a minute and a half now. Maybe I can speed it up.
The executable is called "together".

### Summary Variables
Each pen has some number of susceptile, exposed, infected,
or recovered animals. The coalesced transitions could be faster
if we kept a separate count of these totals per pen. In the
sense that the GSPN is formally a computing language, I can
just make sure that, while the real model computes disease
states, it also keeps current the summaries of the total
numbers in each state per pen.

Turns out this took a bit to get correct, and it's no faster.
This is also in the "together" executable. I overwrote the last
version because I was so sure this would be an excellent improvement.
Thanks to Git for keeping version history.

### Riders
Let's make the first version we hope will scale well.
Instead of pretending each cattlebeast is in mysterious contact with
each other cattlebeast, let's model directly the transfer
of disease. This helps scalability because each animal produces
disease, which then travels, and then infects the other animal.
Instead of having N^2 transitions, the total number of
transitions will be the number of pens times the pen size squared.

This work is in rider.cpp. I made a "pen rider," which I know
isn't always someone on a horse, but you get the point.
"Pen walker" is less fun. The rider walks the sequence of pens
in order, possibly picking up disease and infecting another animal
in this or one of the next pens.

Got this implemented by the end of the day, but now the
simulation doesn't stop, of course. Even if every cattlebeast
is infected, the rider will keep riding infinitely. Have to
tackle that tomorrow. Also, need to come up with some
nice way for the rider to take that ride in the morning
and start on a ride the next day.

## Wednesday 15 October 2014

Late start (11am), and I have a lunch date. The first goal is
to stop the simulation. Then want to put in some reasonable
rates for transitions. Then write the trajectory, meaning
the sequence of events, so that Chris can take a look at it
and make a visualization. Along the way, let's see how
the scaling is, meaning let's see how many animals we can
simulate in 30s, five minutes, or an hour.

Figured out how to get the simulation to stop. Haha.
With 32 pens, 3200 animals not leaving the pen, and a rider
that moves through every pen once a day, for 14000 transitions,
it takes 23s. Used 20% of 16GB memory. Time for lunch.
6400 animals, 32 pens, 25000 transitions, 100s. Used 60% of memory.

- Generate data file for Chris to make visualizations of pens.
- Find why memory usage is high.
- Go back to fused transitions because they will lead to much larger
  herd sizes. Unclear why they seemed slow.

Just spent two hours trying to get the Parameter class to work
correctly, given that some files are shared among applications
and some aren't. Chagned the HDF file to use function templates.
C++ is just a pain all the time.

## Thursday 16 October 2014

### Memory usage is in graph conversion
Thanks to sleeplessness following a 4am fire alarm, I know that
the large memory usage happens during transformation of the
GSPN graph from a list structure to a vector-based structure.
That should, fortunately, be simple to fix with a by-hand rewrite.
Let's first go to scalability, though, by looking again at coalesced
transitions.

There were two ways to check. One was valgrind, which worked
well immediately. The other is gperftools' use of tcmalloc,
which gave me an error about headers >=2**16. Dunno.

### Scalability surprisingly poor for coalesced

Individuals | Pens | Per pen | Expanded [s] | Coalesced [s]
----------- | ---- | ------- | ------------ | -------------
1024        | 32   | 16      | 18           | 10
2048        | 32   | 32      | 81           | 60
4096        | 32   | 64      | tl;dr        | 275
2048        | 64   | 16      | 81 (tl;dr)   | 43
4096        | 256  | 16      | tl;dr        | 180

That looks like the coalesced do worse for large numbers. My
implementation kept separate places for each individual, though.
I could combine those places and put all the susceptible tokens
at the same place, for this model, at least.

Scalability comes from a couple of factors. 1) How many
transitions are affected each time one transition fires?
This determines how many updates happen at each step.
2) How many total transitions are active at any time?
This determines how long it takes to select the next
transition from the list of times. What I'm seeing is that, during
the fastest part of the epidemic, there about sqrt(N)
infecteds and sqrt(N) susceptibles, and, even for the coalesced
transitions, every neighboring transition on the GSPN
has to be checked 

### Put all tokens for S and R on same place

Let's get rid of individual S states and R states, putting all tokens
on one place per pen.

Individuals | Pens | Per pen | Expanded [s]  | Single-SR [s]
----------- | ---- | ------- | ------------  | ---------
1024        | 32   | 32      | 18                        | 5.7
2048        | 32  | 64      |  81                      | 30
4096        | 32  | 128      |  swap                      | 190
2048        | 128  | 16      |  81                      | 21
4096        | 128   | 32      | tl;dr                  | 100
2048        | 256   | 8      |  tl;dr                | 23
4096        | 512  | 8      | tl;dr         | 87
4096        | 256  | 16      | tl;dr        | 91

I realized the problem. When checking which transitions are enabled,
I was looking at the transition outputs. I modified that. Now
it's 8.3s for 2048, 35s for 4096. Of these, the time spent in the
inner loop is 2.8s for 2048 and 12.4s for 4096. That's good.
Time to return to the memory problem.

### Tackling memory usage

Memory usage came from building the GSPN. I made the ability to
bypass the super-friendly gspn-builder in favor of building
the GSPN directly. The penalty you pay is that it could be slow
adding places and transitions while it reallocates memory, but
if you know beforehand about how many places and transitions you need,
it should be faster to build and will, in all cases, use less memory
during the build process.

Using the Single S and R code, timings and transition counts are as
follows. Note that this code creates many more transitions than the
equivalent rider code, which is closer to what we will use.

Individuals | Pens | Places | Transitions | Mem usage [Gb] | Single-SR [s]
----------- | ---- | ------ | ----------- | -------------- | -------------
1024 | 32 | 2176 | 34816 | 0.144 | 1.9
2048 | 32 | 4224 | 69632 | .432 | 7.5
4096 | 32 | 8320 | 139267 | 1.63 | 31
8192 | 32 | 16512 | 278528 | 5.84 | 130

The SingleSR code makes N(P+2) transitions. The rider code makes
2N transitions for the E-I and I-R transitions and approximately another
3P(N/P) transitions for infection of other cattle. That's total=N(2+3/P).
Memory usage comes from both the marking and storage of transitions,
so we'll estimate there should be approximately ten times more cattle we can
fit into memory with the rider graph.

## Friday 17 October 2014

### Transfer improvements into the rider code

Well, that took an hour. The rider now uses coalesced transitions.
Please forgive the powers of two. There is no restriction to powers
of two, just habit.

Individuals | Pens | Per Pen | Transitions | Mem usage [Gb] | Rider-SR [s]
----------- | ---- | ------ | ----------- | -------------- | -------------
1024 | 32 |  32  | 0.032 | 0.90
2048 | 32 |  64  | 0.064 | 3.2
4096 | 32 |  128  | 0.176 | 14.5
8192 | 32 | 256   | 0.592 | 67
16384 | 32 | 512   | 2.032  | 330
16384 | 64 | 256   | 1.216  | 95
16384 | 128 | 128   | 0.816  | 38
32678 | 256 | 128 | 2.160 | 72
65536 | 512 | 128 | 5.824  | 220
102400 | 1024 | 100 | 12  | 560

Now that we have scalability, it's time to get the distributions
correct. They include a Weibull for the latent period and
a Gamma for the infectious period. I've implemented the Weibull
in the Semi-Markov library already, but the Gamma is incomplete
(get it? incomplete? get it?).

## Saturday 18 October 2014

### Test individual transitions

Spent too many hours getting the Gamma math straight, but it
matched the theoretical distribution the first time when
I typed it in, so that felt good. Now I'm building an executable
which will run a single individual over and over again in order
to see that it does follow the EIR part of SEIR that was described.
The executable is in individual.cpp. I run it and make plots
with individual.jl.

I also put the model parameters and EIR part of the transitions
into separate header files in order to share them among
executables. They can certainly be read from a file but are
in a header file for now. Refactoring this took a little bit.

## Sunday 19 October 2014

### Put non-exponential transitions into other models.

First, I noticed that Mardones et al have graphs of their
distributions and that those graphs don't match the plots
of the theoretical distributions. Namely, the latent period
in the paper's plot has a higher, sooner peak. Ugh.
I can't fix that right now. It means mailing the authors.

The nonexponential distributions are in the code now,
but the rider implementation is so dorky that nothing ever
finishes. That is, the rider takes 24 hours to travel the
pens, so it gets infected and recovered over and over in each
pen, and the simulation doesn't proceed. I'll commit
the new distributions and go back and fix the rider.

Somehow I introduced a bug. The well-mixed code,
"together", works fine. The rider code stops
after one transition. Can't figure it out. Going home.

## Monday 20 October 2014

### Solve the rider bug

The list of initially-enabled transitions looks correct.
The times aren't infinite, so something should be enabled.
There is also a last infected that doesn't recover.

Here's how I solved the problem:

- Added a validation test for Uniform distributions to the
  test suite, as I should have done from the start.
  This allowed me to fix sampling uniform distributions.

- Created a check on the invariants for the marking. This
  lead to noticing that one of the transitions failed to
  move a token.

- Made a way to peek at individual transitions associated
  with places. This speaks to a lack of an API with which
  to debug your GSPN.

I also tried lots of other things, including printfs in
places, or walking through with the visual debugger.
These were not helpful in this case.

Those steps lead me to realize that, if we put all S
tokens in the same place, it's hard to know to which E place
to move the S token when we take it, so I've now labeled
all S tokens by their index within the pen.

## Tuesday 21 October 2014

### Format for saved files

I'm realizing that plotting all of this, with the pens,
is complicated, and that the file format I devised is
incomplete. All the time comes down to things like file formats,
doesn't it?

