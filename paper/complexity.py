#!/usr/bin/env python

from rearr import complexity, config

if __name__ == "__main__":
    cfg = config.complexity_config()

    complexity.generate_complexity_ref(cfg)
    complexity.generate_complexity_edit(cfg)
