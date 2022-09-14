#!/usr/bin/env python

import models

engine = models.get_engine()

models.Base.metadata.create_all(engine)

