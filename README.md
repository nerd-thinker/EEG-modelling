# From Brainwaves to Biomarkers
### A feature-extraction pipeline for predicting cognitive age from resting-state EEG

This repository contains the code, references, and acknowledgments for the SPUR poster
*"From Brainwaves to Biomarkers"*.

---

## Feature extraction — core function

Run once per band × channel × scaling combination (5 bands × 32 channels × 4 scalings
per participant):

```r
# Full code can be found in ~/eeg/feature_engineering/complete_feature_function.R

get_features_single <- function(x, time, scaling = "z_score", k = 70) {
  norm_x   <- scale(x)                        # normalise
  fit      <- gam(norm_x ~ s(time, k = k))    # GAM smooth
  smooth_x <- predict(fit, data.frame(time = time))

  peaks <- findpeaks(smooth_x,
    minpeakheight   = mean(smooth_x) + sd(smooth_x),
    minpeakdistance = 10)

  c(mean = mean(smooth_x),
    skewness   = skewness(smooth_x),
    n_peaks    = nrow(peaks),
    activity   = hjorth(smooth_x)$activity,
    mobility   = hjorth(smooth_x)$mobility,
    complexity = hjorth(smooth_x)$complexity,
    sample_entropy   = sample_entropy(smooth_x),
    spectral_entropy = spectral_entropy(smooth_x))
}
```

**k = 70** sets the maximum number of basis functions in the GAM smooth. The fit is then penalised back down
via cross-validation, so k=70 gives enough room to trace real oscillations without
letting the curve chase sample-to-sample noise.

---

## References
Chaudhary, U. (2025). Machine learning with brain data. In U. Chaudhary (Ed.), 
_Expanding senses using neurotechnology: Volume 1 – Foundation of brain-computer 
interface technology_ (pp. 179–223). Springer Nature Switzerland. https://doi.org/10.1007/978-3-031-76081-5_5

Tadel, F., Baillet, S., Mosher, J. C., Pantazis, D., & Leahy, R. M. (2011). Brainstorm: A 
user-friendly application for MEG/EEG analysis. _Computational Intelligence 
and Neuroscience_, 2011, Article 879716. https://doi.org/10.1155/2011/879716

Thornberry, C., Fox, R., Wozniak, A., & Commins, S. (2026). Age-related resting state 
EEG differences in learning and memory performance during a spatial learning task. 
_Brain and Cognition_, 193, Article 106393. https://doi.org/10.1016/j.bandc.2025.106393

---

## Acknowledgements

With thanks to SPUR, the Hamilton Institute for financial support, and the Department
of Mathematics and Statistics for technical support; Department of Psychology for
access to data. With thanks also to Rafael Moral, Sean Commins, and Conor Thornberry
for mentorship and support throughout this project.
