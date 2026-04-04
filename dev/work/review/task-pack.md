---
job: review
mode: persistent-workspace
intent: review
project: ipdcstorm
sources:
  - id: primary_source
    path: sources/source.md
stages:
  - id: draft
    role: technical-writer
    order: 1
    input: source
    goal: Describe the task for this stage.
deliverables:
  - name: final
    from: draft
    path: outputs/final.md
---

# Task

This task pack belongs to a persistent workspace under `dev/work/<intent>`.
Reuse it for ongoing work in the same intent area.

# Objective

Describe the successful outcome.
