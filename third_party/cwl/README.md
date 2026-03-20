# BD Rhapsody CWL pipelines

The files in this directory are **external to rhapsodist** and are not developed
or maintained here. They are the official BD Rhapsody analysis pipeline workflows
published by BD Genomics / CRSwDev at:

  https://bitbucket.org/CRSwDev/cwl

All copyright and licensing rights remain with their authors. Please refer to the
Bitbucket repository for the applicable terms.

We are grateful to the BD Genomics team for publishing the pipeline as open,
human-readable CWL and for keeping it up to date. Following the CWL standard
makes it possible to run the official pipeline alongside open-source aligners in
a reproducible way, which is the whole point of rhapsodist.

## Versions included

| Version | File | Source |
|---------|------|--------|
| 2.2.1 | `v2.2.1/rhapsody_pipeline_2.2.1.cwl` | https://bitbucket.org/CRSwDev/cwl/src/master/v2.2.1/ |
| 3.0   | `v3.0/rhapsody_pipeline_3.0.cwl`     | https://bitbucket.org/CRSwDev/cwl/src/master/v3.0/   |

## Retrieving a version

```bash
curl -sL https://bitbucket.org/CRSwDev/cwl/raw/master/v2.2.1/rhapsody_pipeline_2.2.1.cwl \
     -o docker/cwl/v2.2.1/rhapsody_pipeline_2.2.1.cwl

curl -sL https://bitbucket.org/CRSwDev/cwl/raw/master/v3.0/rhapsody_pipeline_3.0.cwl \
     -o docker/cwl/v3.0/rhapsody_pipeline_3.0.cwl
```

Point `sbg_cwl` in your config to the version you want to use.
