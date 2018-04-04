[![Build Status](https://travis-ci.org/OpenHumans/oh-ancestrydna-source.svg?branch=master)](https://travis-ci.org/OpenHumans/oh-ancestrydna-source)
[![Test Coverage](https://api.codeclimate.com/v1/badges/270ba48837d32c387e6d/test_coverage)](https://codeclimate.com/github/OpenHumans/oh-ancestrydna-source/test_coverage)
[![Maintainability](https://api.codeclimate.com/v1/badges/270ba48837d32c387e6d/maintainability)](https://codeclimate.com/github/OpenHumans/oh-ancestrydna-source/maintainability)

# The *AncestryDNA* upload project for *Open Humans*

**WORK IN PROGRESS**

[![Deploy](https://www.herokucdn.com/deploy/button.svg)](https://heroku.com/deploy)

This is a Django project that is based on the [23andMe Uploader](https://www.github.com/OpenHumans/oh-23andme-source) (which in turns is based on [Open Humans Data Uploader](https://www.github.com/gedankenstuecke/oh_data_uploader) and [FamilyTreeDNA Uploader](https://www.github.com/gedankenstuecke/ftdna-upload)). It uses the same general logic for setting up the project. It should come with a *Celery* task that is enqueued as soon as a new file is uploaded.

This task grabs the newly uploaded file from *Open Humans* and performs the file format verifications to test that it's a proper *23andMe* file. Invalid rows from the file are being deleted to ensure the privacy of participants. If the uploaded file is completely broken the file will be deleted and the user is getting a message notifying them about the problem.

The integration of the *VCF* conversion is yet to be done.
