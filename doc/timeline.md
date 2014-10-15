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

