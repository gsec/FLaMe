Some kind of HEADING
====================

This is the main documentation file. At least so far.

Goals
~~~~~

This is my thesis project. Here we have two main aspects: **Science** and
**Coding**.  While I am studying physics and therefore the science aspect is
important, equally, if not even more valuable, is to be able to produce good
quality software.

We will be using `Python` as language due to it's expresiveness and simplicity.
It is very powerful supporting different programming paradigms and the ability
to express code in a very clear, compact and human-readable style.

Scientific Aspects
------------------

Scientific question.
  What exactly are we trying to to and why? Here we discuss the purpose of the
  simulation and the relevance to other projects and goals we have in our
  science group.

Simulation development.
  This will be a journal of the development and especially on the refinements of
  the methods and on the thinking process of our simulation.  As we start, we
  have a very big and vague picture in our mind how the simulation will look
  like and which features can be thought of. Therefore we will start with a very
  basic and rough framework. This framework should be simple but general. This
  enables us to refine and improve the model we are working with. Up to a degree
  we think is useful.  This approach will enable a freedom to go in whichever
  direction we might want, as ideas as well as problems come along the way and
  can not always be foreseen.

Related fields.
  At the beginning of the work in our group we examined various physical systems
  closely related to the flakes we want to reproduce here. There should be some
  space to address those topics, depending on the relevance and potential
  advancements thereof.  Some examples might be:

  * Plasmons
  * Nano wires
  * Quantum dots
  * Signal enhancement
  * Amplification and lasing
  * Chemical and physical growth aspects
  * Environment of the flake and pre-conditions of the growth

Coding Aspects
--------------

Clear code.
  We want the code to be self-explanatory as far as possible. The instructions
  should be as close to the thinking process as possible.

No repetitions.
  The same code written many times in slight variations only introduces typos
  and mistakes. It makes the code hard to maintain and hard to follow. Therefore
  we will avoid repetitions and use `refactoring` as a tool to maintain sleek
  code as we advance with our project.

Testing & Coverage.
  Tests are essential in software production. Code that is not tested can not be
  relied on, independent of how well-written it actually is. Code which has no
  tests, can not be improved in safe manners guaranteeing unaltered
  functionality. We strive here for full coverage with the usage of unit tests.

Documentation.
  This is important on different layers. At first for the main-developer
  himself.  When you see code that has been written as short as a few months
  ago, one can already struggle in understanding certain functions, variables
  and reasoning for certain ideas.  Only because ideas can not be maintained
  fresh for a long time and present for all times, as soon as they get slightly
  complex.

  In the next level documentation is important for collaborators. It prevents
  misuse and correct and efficient usage of the software. This applies even more
  in a environment of open source software where projects and collaborators can
  be all over the world and where a face to face communication is not an option.

  At last documentation is indispensable for the end-user. Often the
  documentation of a software will be the first, and sometimes the only, source
  of informations about how to use and the purposes of a program.

  On the one hand we have comments and doc strings in our source code itself. On
  the other hand we are building a documentation independent from that by using
  `Sphinx` to give a deeper insight in the ideas and workings of the project.
  This will also be used to document the process of thought, and if feasible,
  also be used to write the thesis.
