Goals
~~~~~

This is my thesis project. Here we have two main aspects: **Science** and
**Coding**. While I am studying physics and therefore the science aspect is
important, we also aim to understand the tools and methods behind proper
software development.

We will be using `Python` as language due to it's expressiveness and simplicity.
It is a powerful language, supporting different programming paradigms and the
ability to compose code in a very clear, compact and human-readable style.

Scientific Aspects
------------------

Scientific question.
  What exactly are we trying to realize and why? Here we discuss the purpose of
  the simulation and the relevance to other projects and goals we have in our
  science group.

Simulation development.
  This will be a journal of the development and especially on the refinements of
  the methods and on the thinking process of our simulation.  As we start, we
  have a very big and vague picture in our mind how the simulation will look
  like and which features can be thought of. Therefore we will start with a very
  basic and rough framework. The framework should be simple but general.
  Enabling us to refine and improve the model we are working with, up to a
  degree we think is useful. This approach will enable a freedom to go in
  whichever direction we might want. Ideas, as well as problems, come along the
  way and can not always be foreseen.

Related fields.
  At the beginning of the work in our group we examined various physical systems
  closely related to the flakes we want to reproduce here. There should be some
  space to address those topics, depending on the relevance and potential
  advancements thereof. Some examples might be:

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
  We want the code to be self-explanatory as far as possible. The functions
  and instructions should follow the thought process and be kept short.

No repetitions.
  The same code written many times in slight variations introduces typos and
  mistakes. It makes the code hard to maintain and hard to follow. Therefore we
  will avoid repetitions and use `refactoring` as a tool to maintain sleek code
  as we advance with our project.

Testing & Coverage.
  Tests are essential in software production. Code that is not tested can not be
  relied on, independent of how well-written it actually is. Code which has no
  tests, can not be improved in safe manners guaranteeing unaltered
  functionality.

Documentation.
  This is important on different layers.
  At first for the main-developer himself. When you look at code written a few
  months ago, one already struggles in understanding certain functions,
  variables and reasoning for certain ideas.

  In the next level, documentation is important for collaborators. It is
  essential for correct and efficient usage of the software. Ideas shared at
  this level ideally reveal a broad picture of the program. This applies even
  more in a environment of open source software. Projects and collaborators can
  be all over the world and face to face communication not always an option.

  At last documentation is indispensable for the end-user. Often the
  documentation of a software will be the first, and sometimes the only, source
  of informations about usage and purpose the program.

  On the one hand we have comments and doc strings in our source code itself. On
  the other hand we are building a documentation independent from that by using
  `Sphinx` to give a deeper insight in the ideas and workings of the project.
  This will also be used to document the process of thought, and if feasible,
  also be used to write the thesis.
