Contributing
============

Thank you for your interest in contributing to SpectralUnmxing. If you are just getting
started, please review the guidelines below to understand how and where you can
best support and contribute to this project.  Typical contributions may include:

* Code patches
* Feature enhancements
* Documentation improvements
* Bug reports

If you have ideas for new additions, that's great - please contact the maintainers
at the addresses given below, and we can coordinate efforts.  Our general policy
is to for the maintainers to delegate technical authority to individuals to make
changes and additions to specific topics.


Getting Started
---------------

First of all, to get a sense of the project's current status and roadmap, please
be sure to spend some time reviewing issues in the `issue tracker <https://github.com/emit-sds/SpectralUnmixing/issues>`_.

If you have have discovered a new issue or task, then go ahead and create a `new
issue <https://github.com/emit-sds/SpectralUnmixing/issues/new>`_.


Fork and Create a Branch
------------------------

SpectralUnmixing follows the `Standard Fork and Pull Request <https://gist.github.com/Chaser324/ce0505fbed06b947d962>`_ workflow.

When you have chosen an issue to work on, start by `forking <https://help.github.com/articles/fork-a-repo/>`_ the repo.

Then, create a branch with a descriptive name.  A good branch name starts with
the issue number you are working on followed by some descriptive text.  For
example, if you are working on issue #113 you would run:

```
  git checkout -b 113-add-spectral-residual
```


Testing
-------

Tests live in `test/ <test/>`_, and are executed using
```
julia --project=@. -e 'using Pkg; Pkg.test()'
```

Our development strategy employs continuous integration and unit testing to validate all changes.  We appreciate you writing additional tests for new modifications or features.  Depending on the significance of your changes, additional tests on real data may be requested.


Implement Your Changes and Create a Pull Request
------------------------------------------------

At this point, you are ready to implement your changes!

As you develop, you should make sure that your branch doesn't veer too far from SpectralUnmixing's dev branch.  To do this, switch back to your dev branch and make
sure it's up to date with SpectralUnmixing's dev branch:

```
  git remote add upstream https://github.com/emit-sds/SpectralUnmixing.git
  git checkout dev
  git pull upstream dev
```

Then update your feature branch from your local copy of dev, and push it!  We recommend using git's rebase call when possible.

```
  git checkout 113-add-spectral-residual
  git rebase dev
  git push --set-upstream origin 113-add-spectral-residual
```

When you are ready to submit your changes back to the SpectralUnmixing repo, go to GitHub
and make a `Pull Request <https://help.github.com/articles/creating-a-pull-request/>`_

Keeping your Pull Request Updated
---------------------------------

If a maintainer asks you to "rebase" your PR, they're saying that a lot of code
has changed, and that you need to update your branch so it's easier to merge.

Here's the suggested workflow:

```
  git checkout 113-add-spectral-residual
  git pull --rebase upstream dev
  git push --force-with-lease 113-add-spectral-residual
```

Project Decision Making
-----------------------

Minor changes follow an expedited acceptance process.  These are things like:

* Bug fixes
* Unit tests
* Documentation
* Consolidation that does not change algorithm results or provide significant new functionality
* New functionality initiated by maintainers, or over which authority has been delegated in advance by maintainers (e.g. through issue assignment)

Minor change pull requests are accepted once they:

* Pass unit tests and adhere to project coding conventions
* Get signoff from at least one maintainer, with no objections from any other maintainer

Accepted minor changes will be released in the next major or minor release version. Hotfixes will be expedited as needed.

Major changes include:

* New functionality, including examples, data, and algorithm changes, over which authority was not delegated in advance.
* Official releases
* Project policy updates

These are accepted through consensus of a quorum of maintainers.  **If you would like to include any new algorithms or examples, we highly recommend that they are supported by peer reviewed scientific research.**

Release Steps (for Maintainers)
-------------------------------

Releases should trigger a new PyPi upload, and subsequently a fresh upload to conda-forge.  Therefore,
the revised steps for versioning are:

* Submit version number change to setup.cfg in dev
* Trigger a PR from dev to main
* Accept the PR
* Go to https://github.com/isofit/isofit/releases
* Click "Draft a new release"
* Enter tag version as "v3.8.0" (depending on latest version), and input release title and description
* Click "Publish release"

Contributors
------------

The github maintainers, responsible for handling pull requests, are:

* Philip G. Brodrick philip.brodrick@jpl.nasa.gov
* Vatsal Jhalani (NASA JPL)

Thanks to the following regular contributors:

* Francisco Ochoa (UCLA)

The SpectralUnmixing codebase was made possible with support from the Earth Surface Mineral Dust Source Investigation (EMIT), an Earth Ventures-Instrument (EV4) Mission.
The initial research took place at the Jet Propulsion Laboratory, California Institute of Technology, 4800 Oak Grove Dr., Pasadena, CA 91109 USA.
