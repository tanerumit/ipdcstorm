# Tutorial 3 - Climate Impact Analysis

> **Prerequisites**
>
> This tutorial assumes you have already worked through [Tutorial
> 1](https://tanerumit.github.io/ipdcstorm/articles/tutorial_1_introduction_setup.qmd),
> where the single-site workflow and the main output objects are
> introduced in more detail.

This tutorial shows how to move from one island to a small group of
islands, then compare a stationary baseline with a warmer-climate
scenario. The focus is practical: run the model, inspect the annual
outputs, and then downscale one site to daily disruption metrics.

## 1 What This Tutorial Covers

By the end, you will be able to:

- run the baseline model for three islands at once
- compare a stationary baseline with a warmer climate scenario
- translate annual scenario differences into daily disruption metrics
  for one site

``` r
library(ipdcstorm)
library(dplyr)
library(ggplot2)
```

## 2 Multi-Site Baseline

We start with the three Dutch Caribbean islands used throughout this
package: St. Martin, Saba, and Statia. The setup is the same as in
Tutorial 1, except that `targets` now contains three rows instead of
one.

``` r
seed <- 123L
simulation_years <- 500L

ibtracs_demo <- system.file(
  "extdata",
  "ibtracs_demo.csv",
  package = "ipdcstorm"
)

targets <- data.frame(
  name = c("St_Martin", "Saba", "Statia"),
  lat = c(18.0708, 17.6350, 17.4890),
  lon = c(-63.0501, -63.2300, -62.9740)
)

cfg_baseline <- make_hazard_cfg(
  data_path = ibtracs_demo,
  historical_start_year = 1970L,
  search_radius_km = 800,
  simulation_years = simulation_years,
  climate = make_climate_cfg(scenario = "stationary")
)
```

``` r
out_baseline <- run_hazard_model(
  cfg = cfg_baseline,
  targets = targets,
  seed = seed,
  verbose = FALSE
)
```

The `rates` table is the historical calibration result: mean annual
site-level rates for tropical storms and hurricanes at each island.

``` r
out_baseline$rates
```

    #> # A tibble: 6 × 6
    #>   location  storm_class lambda n_years prob_annual prob_none
    #>   <chr>     <chr>        <dbl>   <int>       <dbl>     <dbl>
    #> 1 St_Martin HUR          0.191      47       0.174     0.826
    #> 2 St_Martin TS           0.553      47       0.425     0.575
    #> 3 Saba      HUR          0.128      47       0.120     0.880
    #> 4 Saba      TS           0.511      47       0.400     0.600
    #> 5 Statia    HUR          0.149      47       0.138     0.862
    #> 6 Statia    TS           0.511      47       0.400     0.600

The annual simulation table gives one synthetic season per row. For a
quick comparison across islands, we can average the simulated counts and
hurricane share by location.

``` r
baseline_summary <- out_baseline$sim |>
  group_by(location) |>
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_ts = mean(n_ts, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    mean_p_hur = mean(p_hurricane, na.rm = TRUE),
    .groups = "drop"
  )

baseline_summary
```

    #> # A tibble: 3 × 5
    #>   location  mean_total mean_ts mean_hur mean_p_hur
    #>   <chr>          <dbl>   <dbl>    <dbl>      <dbl>
    #> 1 Saba            1.19   1.02     0.17       0.141
    #> 2 St_Martin       1.01   0.814    0.194      0.186
    #> 3 Statia          1.25   1.09     0.16       0.137

| Column | Description |
|----|----|
| `mean_total` | Mean annual number of site-impacting storms with at least tropical-storm strength. |
| `mean_ts` | Mean annual number of tropical-storm events. |
| `mean_hur` | Mean annual number of hurricane events. |
| `mean_p_hur` | Mean simulated fraction of storms that are hurricanes. |

These averages should broadly line up with the rate table. Small
differences across islands reflect real differences in how often storms
pass close enough, and strongly enough, to exceed the site thresholds.

## 3 Climate Change Scenarios

The climate extension changes the annual simulation in two SST-driven
ways. First, warming can scale the total number of storms up or down
through a rate effect. Second, warming can shift the storm mix toward a
larger hurricane fraction through an intensity effect.

In this model, both effects are driven by a scalar `delta_sst`: the Main
Development Region sea-surface temperature change, in °C, relative to a
1991-2020 baseline. That keeps the workflow simple for scenario testing:
you can either use a named scenario helper or supply a direct warming
value yourself.

We will compare a stationary baseline (`delta_sst = 0`) with a
hypothetical warmer climate using `delta_sst = 1.5`. That is a clean
side-by-side example for understanding the direction and size of the
climate response.

``` r
sst_scenario_info() |>
  filter(source == "ipcc_ar6") |>
  select(scenario, source, delta_sst_2050, delta_sst_2100)
```

    #> # A tibble: 3 × 4
    #>   scenario source   delta_sst_2050 delta_sst_2100
    #>   <chr>    <chr>             <dbl>          <dbl>
    #> 1 ssp126   ipcc_ar6            0.3            0.4
    #> 2 ssp245   ipcc_ar6            0.5            1  
    #> 3 ssp585   ipcc_ar6            1              2.5

The built-in scenario table is useful for context, but for this tutorial
we set the warming directly so the climate signal is explicit in the
code.

``` r
climate_future <- make_climate_cfg(
  delta_sst = 1.5,
  sensitivity_mode = "fixed"
)

cfg_future <- make_hazard_cfg(
  data_path = ibtracs_demo,
  historical_start_year = 1970L,
  search_radius_km = 800,
  simulation_years = simulation_years,
  climate = climate_future
)
```

``` r
out_future <- run_hazard_model(
  cfg = cfg_future,
  targets = targets,
  seed = seed,
  verbose = FALSE
)
```

Now we stack the two annual simulation tables and compute a compact
comparison by scenario and island. The summary below reports the mean
annual storm count, the tropical-storm and hurricane components, and the
mean simulated hurricane fraction.

``` r
activity_summary <- bind_rows(
  out_baseline$sim |> mutate(scenario = "Baseline"),
  out_future$sim |> mutate(scenario = "Future (+1.5 C)")
) |>
  group_by(scenario, location) |>
  summarise(
    mean_total = mean(n_total, na.rm = TRUE),
    mean_ts = mean(n_ts, na.rm = TRUE),
    mean_hur = mean(n_hur, na.rm = TRUE),
    mean_p_hur = mean(p_hurricane, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(location, scenario)

activity_summary
```

    #> # A tibble: 6 × 6
    #>   scenario        location  mean_total mean_ts mean_hur mean_p_hur
    #>   <chr>           <chr>          <dbl>   <dbl>    <dbl>      <dbl>
    #> 1 Baseline        Saba            1.19   1.02     0.17       0.141
    #> 2 Future (+1.5 C) Saba            1.24   1.03     0.204      0.168
    #> 3 Baseline        St_Martin       1.01   0.814    0.194      0.186
    #> 4 Future (+1.5 C) St_Martin       1.15   0.91     0.238      0.220
    #> 5 Baseline        Statia          1.25   1.09     0.16       0.137
    #> 6 Future (+1.5 C) Statia          1.31   1.11     0.202      0.162

| Column       | Description                                               |
|--------------|-----------------------------------------------------------|
| `scenario`   | Baseline or future climate label used for the comparison. |
| `location`   | Target island.                                            |
| `mean_total` | Mean annual storm count under that scenario.              |
| `mean_ts`    | Mean annual tropical-storm count.                         |
| `mean_hur`   | Mean annual hurricane count.                              |
| `mean_p_hur` | Mean simulated hurricane fraction.                        |

``` r
activity_summary |>
  ggplot(aes(x = location, y = mean_total, fill = scenario)) +
  geom_col(position = "dodge") +
  labs(
    x = NULL,
    y = "Mean annual storm count",
    fill = "Scenario",
    title = "Annual storm activity by scenario and location"
  ) +
  theme_minimal(base_size = 11)
```

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABUAAAAPACAMAAADDuCPrAAABoVBMVEUAAAAAADoAAGYAOjoAOmYAOpAAZrYAv8Q6AAA6OgA6Ojo6OmY6OpA6ZmY6ZpA6ZrY6kJA6kLY6kNtNTU1NTW5NTY5Nbm5Nbo5NbqtNjshmAABmADpmOgBmOjpmOpBmZgBmZjpmZmZmZpBmkGZmkJBmkLZmkNtmtrZmtttmtv9uTU1uTY5ubk1ubm5ubo5ujo5ujqtujshuq6tuq+SOTU2Obk2Obm6Ojm6Ojo6Oq6uOq8iOyOSOyP+QOgCQOjqQZjqQZmaQZpCQkGaQkLaQtraQttuQ29uQ2/+rbk2rjm6rq46rq8iryOSr5Mir5P+2ZgC2Zjq2Zma2kDq2kGa2kJC2tma2tpC2tra2ttu229u22/+2///Ijk3Ijm7Iq27Iq47Iq6vIyKvI5KvI5OTI5P/I/8jI///bkDrbkGbbtmbbtpDbtrbbttvb27bb29vb2//b/7bb///kq27kq47kyI7kyKvkyMjk5Kvk5Mjk5P/k/+Tk///r6+v4dm3/tmb/yI7/25D/27b/29v/5Kv/5Mj/5OT//7b//8j//9v//+T///8mv9pzAAAACXBIWXMAAB2HAAAdhwGP5fFlAAAgAElEQVR4nO3di58c55Xe9xoRQIagQXoxlETR4AaQyTWBrNdSRNmJLYSMdmlRvuxwSccSKXptJ8yOAC1Fg5Hii5jBhTQ4mL86XdWnu6vqfavnqe7TdWaqft/PZ1fAYC41T53zsHu6u6Y4BQBspIg+AAC4qChQANgQBQoAG6JAAWBDFCgAbIgCBYANUaAAsCEKFAA2RIECwIYoUADYEAUKABuiQAFgQxQoAGyIAgWADVGgALChwQr0sCj2fraTz/z01V19Zh9dx9fnuE/uFsVtz4PysIvgZ2NyfRdf5LwPCS6ooQp0NsDt1XD81MJu3H9t8P2xL3l2gZ59bBTops44CcBWhirQ41l/Fpcf7OJTK7sxK6Ch92f5Jc8qUOXYKNDNnHkSgK0MVaCHxd7zO+oAZTcC9kf+kuJ/AM5fge6Cc4HSm9itgQr00X5x+V5RXNnF56ZAR4QCxYUyUIEezfZiV8NMgY4IBYoLZZgCnf8oKlkOHxToiFCguFCGKdDj6gGk2f//1q8WbzosG+Hx2/tF8dwbD9a9abFPR6ufANwv32f2Tp9Xf8styeO3ny/f48WPqr89Kt+9WD6P6rPXZn/fe/Hj1Se+ffL+frH3Qtnxs/epDuGFd6uv9Nrsw174KPMtNY+hdPJh+TXnH1f/kvPjO6r/AKP6tuZvX71j4/s4bBbmvEA/e6087o/S7/owqddmAOnx2Ts1gsjkn36jtbRqh9CKdF1S2a9yUga999KDtQWafpH0W2p8seQk5D5N/rsGJMMU6HwtGjejyrl9v2gUW/5NSYE+vlEs3qn6bGmBniw+S1FcKvekUaDlUcy99GDxiW8f2j+XBXpk/3x99XmSG87tY5j5tX2R+XMNkt2d/f/lfz2Wb2kUaP2brb/3qR307cWBvfigFU37vZMA0uNrvJMFkck//UZraS2DTyJdl9S6r7L3s+4CzXyR5FtqfbFcgSafJnc8gGiQAl1M71HtmUyzub263PLLD7rf1C7Q6hml9YFPC/Rw9R5Vs9QLtP7h869xtPiyl8sbQHs/Xf7zm8vP016s5Biqz1L/tOnu1m4lzm+QJwVaPtL2YPkejQfcyrW/1TzuNe+dBJAeX71JFh+dyT/9RmtpLYJPI12XVOarLDKYHez/2FWgmS+SfEvtL5Yp0PTTZI4HUA1SoEernV9WUbXj5f2o6lbD7e43tQu0fEnTmw/Ku2/7839MCrTcm5fK+3BPPijs41fvc2hf47Mbi3+r9vCPHpw++bx2CPfLu3TVFzr5oEiePpAcQ/U1L71rN+uaX9L+UGu5eZfWdnp548gOMvmR57ztLi+yub7+vTMBpMdXftvlOz1ZviWff+sbraW1OPA00nVJpV+l+ubKQ7lf3YDMF2j6RdJvKf1i7ZOQ+TSZ7xpQDVGgywUv/7BYj8PlqpQ3Cq53v6lVoOU/2ZTPb8ilBVr7cWPtZuv8fR4tVuv05NBunR0VuaM6Luo/R2jeMkmPofY+i/9atHd3VXh2jzsp0Oq5Cqd2lK375HeXLV4/7o73zgSQHF/ze6g+Ps0/843W0rIDz0R6RlKtr1J+gjdX32a2QDNf5IxvKXsS8p+mPXWAaogCXS147WGkw8Yfr3S/qVWgx6s2mw18V4Feaqxx7X0OGx9e7dtR0bhdbIdQ26bjdoFmj2Gxu0kv1n6Acb3+naQFurxXftS+zVs2yyKbxdfqfu80gPT46h9kP11I80+/0Xpaq59NtCNdk1TmqxyuDmVVcLUjz3+R9FvKfLH2Scgca2bqANUQBbqay9rQH2ZvJWXelHkU/nT52eq35Jaql41e+4vPm++6fASp/jjWldPkJ7PpoSa38NJjSDo2U6CLwlscQ1qgixup6ZOW6rfdF6F0v3caQHp89Qfuj+efMJN/8o3W05ofeC7SNUmlX6X9CXIFmvkimcjTL9Y6CbljXftdA+sNUKD1n9uvfk6f/+mm8Ka5k9/+23f2i3yBLh8eeeHPFxuWltXpcgVbt8Wa909LXQVaO4bkYZxMgS4Kb1GkmWM6Wt6nbXVDY/MXB7z2vZsBJMdXfwipKJZVcmbY9X+pP5mgGemapNKv0vgER/kCzXyRTOTpF2udhNyxrv1vNLDeAAV63FzV7icoqQX6WfUcx0q2QOvtcGn+HEH3Am0dQ2b1kgJd3NBbdETmmOyeZ9IizQuOLKqj873TALI/E2j/Z00Ju3+Btj9BrkBX8R7LBZptu/YXo0CxW7sv0Paqrp4zs1mBrp7r11mg5VPEl+/SVVbbFGhyDFKBzgsvuUBQ/ZgO7Y5m+xtqvOlolWH+vdMANi3Q5BvtW6DpJ9hhgaZfjALFbu2+QOvPXapmdL4tmxao/UDg6rXXf/mfOn4GOvfk37xWPQ2wUVLtu4tXTjcp0PQYpAKdF97yrmfumKr+SO+T5+/Cd753GkC2QNtP2VHCPqNA218n8wm2vws//7Dsj2nbX2xdgR6d8YML4Gy7L9Cjxo2SR7WnA/Yr0NVDPtVz/U67H0RaOVk8D/KMB5F6Fmh6DPVbXq0H2ZuPsl9fPXqTK9DqRmrHffL2g0jd750GkB5f8mhNLv/0G00LdO2DSJlPkHwVhweRlpXa/mLag0gUKDa18wJtrP78r4tns2gFemX1ga15P87fhc8t5PJ9jtLnsfQv0OwxLL5m+0H25qPsl//T8g5r9lbx7EbqP8/cJ1/Gdtq4SZ9/70wA6fEdNfvn+mku//QbTQs0F+mapPI1Xf8E2acxpV8k/ZYyX6x9EjLHSoFiCzsv0OPWCyEXf9cKdHmPv3oC95VGHz/qehS+tiWLd1++T/6J9P0KNHMMtScYdu3u/MjKV4q2vkL7B3N/99X0PnntGeb1Lu1470wAyfHVnnG5aOTsbcP2N5oW6Jon0nck1TrLzU+wwRPp599S7ou1T0L+01Cg2NTOC7R9p678D3/HPafMm8pH8C/988VjIvZx1f206hWI+QKdv8av3K3yBXvVlqyWK30t30a3QFvHsHpd4fLFk8sv2fxMyfPQm7fR5w/uJPfJ529evtzx9vr3zgSQHl8ZxKWP7HWPHflnvtG0QDORrksqf+JXnyBfoJkvknxLmS+WnIT001Cg2MKuCzS5U7d45YdWoLUHi1/fX95cWSmXIv0Z6FH9Xa7XPs/t5rNS0xfyaAWaHkPjazYu1nG79Sh7+3no9WNbfJ70WafNi4msKib/3pkAOo+vWJZNx23D1jeaKdA00nVJrT3Leze7CjTzRdrfUuaLJSch/TQUKLaw6wI9Shb8eL7T4pM+ny4udHb70f7yXr296R/cbVyUo2Z5nbPZTZLVV7V788vnuqwuZ9f3aUzJMczca19bbfElW4+yL7/H+vNDV29+tJ/b4sbl7GoNk3/vXADp8Z0sn+lkl7zL/wil9Y1mCjSNdF1SubO8aNC92x2Pwue/SPtbypyW5CQkn4YCxRZ2XKCLh35ab7r8QC3Q05N7V+eXv10UqF0G+IV3F08FzD0K/+TD6hJlL7y7XOh75TOs/6j64/3y6T3P1S+o3PuJ9O1jqI60fXVf+5KtR9lbVwVtHVvHtefnb62uOty4anHnlerTADJXH/6suvrw8i25/JNvNFegSaTrksoX1mfVlas/7noaU8cXaX9LmdOSnoTWp6FAsYWhfisnRE/bF0d2fG8AvijQc6bzJd4O7w3AFwV6vuRfmOnz3gCcUaDnSnl1df3XSvR7bwDeKNDz49Hqtfvu7w1gByjQ82P+HMXuF7Zv894AdoACPT/Ke+TPvXn2+23y3gB2gAIFgA1RoACwIQoUADZEgQLAhihQANgQBQoAG6JAAWBDFCgAbIgCBYANUaAAsCEKFAA2RIECwIYoUADYEAUKABuiQBNffhl9BBcHWenIaowo0ASDriMrHVmNEQWaYNB1ZKUjqzGiQBMMuo6sdGQ1RhRogkHXkZWOrMaIAk0w6Dqy0pHVGFGgCQZdR1Y6shojCjTBoOvISkdWY0SBJhh0HVnpyGqMKNAEg64jKx1ZjREFmmDQdWSlI6sxokATDLqOrHRkNUYUaIJB15GVjqzGiAJNMOg6stKR1RhRoAkGXUdWOrIaIwo0waDryEpHVmNEgSYYdB1Z6chqjCjQBIOuIysdWY0RBZpg0HVkpSOrMaJAEwy6jqx0ZDVGFGiCQdeRlY6sxogCTTDoOrLSkdUYUaAJBl1HVjqyGiMKNMGg68hKR1ZjRIEmGHQdWenIaowo0ASDriMrHVmNEQWaYNB1ZKUjqzGiQBMMuo6sdGQ1RhRogkHXkZWOrMaIAk0w6Dqy0pHVGFGgCQZdR1Y6shojCjTBoOvISkdWY0SBJhh0HVnpyGqMKNAEg64jKx1ZjREFmmDQdWSlI6sxokATDLqOrHRkNUYUaIJB15GVjqzGiAJNMOg6stKR1RhRoAkGXUdWui//h/MhOodxoUATlIKOrHQU6BhRoAlKQUdWOgp0jCjQBKWgIysdBTpGFGiCUtCRlY4CHSMKNEEp6MhKR4GOEQWaoBR0ZKWjQMeIAk1QCjqy0lGgY0SBJigFHVnpKNAxokATlIKOrHQU6BhRoAlKQUdWOgp0jCjQBKWgIysdBTpGFGiCUtCRlY4CHSMKNEEp6MhKR4GOEQWaoBR0ZKWjQMeIAk1QCjqy0lGgY0SBJigFHVnpKNAxokATlIKOrHQU6BhRoAlKQUdWOgp0jCjQBKWgIysdBTpGFGiCUtCRlY4CHSMKNEEp6MhKR4GOEQWaoBR0ZKWjQMeIAk1QCjqy0lGgY0SBJigFHVnpKNAxokATlIKOrHQU6BhRoAlKQUdWOgp0jCjQBKWgIysdBTpGFGiCUtCRlY4CHSMKNEEp6MhKR4GOEQWaoBR0ZKWjQMeIAk1QCjqy0lGgY0SBJigFHVnpKNAxokATlIKOrHQU6BhRoAlKQUdWOgp0jHZdoF8CqEQ3p+l1zDuuh4uPW6AJpkZHVrrzUqDROYwLBZqgFHQXI6v/73yIbk4TfTbGhQJNXIxSOB8uRlbRzWmim9NEn41xoUATF6MUzoeLkVV0c5ro5jTRZ2NcKNDExSiF8+FiZBXdnCa6OU302RgXCjRxMUrhfLgYWUU3p4luThN9NsaFAk1cjFI4Hy5GVtHNaaKb00SfjXGhQBMXoxTOh4uRVXRzmujmNNFnY1wo0ATP19NRoD1ED5SJPhvjQoEmKFAdBdpD9ECZ6LMxLhRoggLVUaA9RA+UiT4b40KBJihQHQXaQ/RAmeizMS4UaIIC1VGgPUQPlIk+G+NCgSYoUB0F2kP0QJnoszEuFGiCAtVRoD1ED5SJPhvjQoEmKFAdBdpD9ECZ6LMxLhRoggLVUaA9RA+UiT4b40KBJihQHQXaQ/RAmeizMS4UaIIC1VGgPUQPlIk+G+NCgSYoUB0F2kP0QJnoszEuFGiCAtVRoD1ED5SJPhvjQoEmKFAdBdpD9ECZ6LMxLhRoggLVUaA9RA+UiT4b40KBJihQHQXaQ/RAmeizMS4UaIIC1VGgPUQPlIk+G+NCgSYoUB0F2kP0QJnoszEuFGiCAtVRoD1ED5SJPhvjQoEmKFAdBdpD9ECZ6LMxLhRoggLVUaA9RA+UiT4b40KBJihQHQXaQ/RAmeizMS4UaIIC1VGgPUQPlIk+G+NCgSYoUB0F2kP0QJnoszEuFGiCAtVRoD1ED5SJPhvjQoEmKFAdBdpD9ECZ6LMxLhRoggLVUaA9RA+UiT4b40KBJihQHQXaQ/RAmeizMS4UaIIC1VGgPUQPlIk+G+Nyrgo0esJN9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhoEB7iB4oE302xoUCTUVPuIk+GwoKtIfogTLRZ2NcKNBU9ISb6LOhuBi/Ajp6oEx0SGaYyZgKCjQVPeEm+mwoKNAeokMyw0zGVFCgqegJN9FnQ0GB9hAdkhlmMqaCAk1FT7iJPhsKCrSH6JDMMJMxFRRoKnrCTfTZUFCgPUSHZIaZjKmgQFPRE26iz4aCAu0hOiQzzGRMBQWaip5wE302FBRoD9EhmWEmYyoo0FT0hJvos6GgQHuIDskMMxlTQYGmoifcRJ8NBQXaQ3RIZpjJmAoKNBU94Sb6bCgo0B6iQzLDTMZUUKCp6Ak30WdDQYH2EB2SGWYypoICTUVPuIk+GwoKtIfokMwwkzEVFGgqesJN9NlQUKA9RIdkhpmMqaBAU9ETbqLPhoIC7SE6JDPMZEwFBZqKnnATfTYUFGgP0SGZYSZjKijQVPSEm+izoaBAe4gOyQwzGVNBgaaiJ9xEnw0FBdpDdEhmmMmYCgo0FT3hJvpsKCjQHqJDMsNMxlRQoKnoCTfRZ0NBgfYQHZIZZjKmggJNRU+4iT4bCgq0h+iQzDCTMRUUaCp6wk302VBQoD1Eh2SGmYypoEBT0RNuos+GggLtITokM8xkTAUFmoqecBN9NhQUaA/RIZlhJmMqKNBU9ISb6LOhoEB7iA7JDDMZU0GBpqIn3ESfDQUF2kN0SGaYyZgKCjQVPeEm+mwoKNAeokMyw0zGVFCgqegJN9FnQ0GB9hAdkhlmMqaCAk1FT7iJPhsKCrSH6JDMMJMxFRRoKnrCTfTZUFCgPUSHZIaZjKmgQFPRE26iz4aCAu0hOiQzzGRMBQWaip5wE302FBRoD9EhmWEmYyoo0FT0hJvos6GgQHuIDskMMxlTQYGmoifcRJ8NBQXaQ3RIZpjJmAoKNBU94Sb6bCgo0B6iQzLDTMZUUKCp6Ak30WdDQYH2EB2SGWYypoICTUVPuIk+GwoKtIfokMwwkzEVFGgqesJN9NlQUKA9RIdkhpmMqaBAU9ETbqLPhoIC7SE6JDPMZEwFBZqKnnATfTYUFGgP0SGZYSZjKijQVPSEm+izoaBAe4gOyQwzGVNBgaaiJ9xEnw0FBdpDdEhmmMmYCgo0FT3hJvpsKCjQHqJDMsNMxlRQoKnoCTfRZ0NBgfYQHZIZZjKmggJNRU+4iT4bCgq0h+iQzDCTMRUUaCp6wk302VBQoD1Eh2SGmYypoEBT0RNuos+GggLtITokM8xkTAUFmoqecBN9NhQUaA/RIZlhJmMqKNBU9ISb6LOhoEB7iA7JDDMZU0GBpqIn3ESfDQUF2kN0SGaYyZgKCjQVPeEm+mwoKNAeokMyw0zGVFCgqegJN9FnQ0GB9hAdkhlmMqaCAk1FT7iJPhsKCrSH6JDMMJMxFRRoKnrCTfTZUFCgPUSHZIaZjKmgQFPRE26iz4aCAu0hOiQzzGRMBQWaip5wE302FBRoD9EhmWEmYyoo0FT0hJvos6GgQHuIDskMMxlTQYGmoifcRJ8NBQXaQ3RIZpjJmAoKNBU94Sb6bCgo0B6iQzLDTMZUUKCp6Ak30WdDQYH2EB2SGWYypoICTUVPuIk+GwoKtIfokMwwkzEVFGgqesJN9NlQUKA9RIdkhpmMqaBAU9ETbqLPhoIC7SE6JDPMZEwFBZqKnnATfTYUFGgP0SGZYSZjKijQVPSEm+izoaBAe4gOyQwzGVNBgaaiJ9xEnw0FBdpDdEhmmMmYCgo0FT3hJvpsKCjQHqJDMsNMxlRQoKnoCTfRZ0NBgfYQHZIZZjKmggJNRU+4iT4bCgq0h+iQzDCTMRUUaCp6wk302VBQoD1Eh2SGmYypoEBT0RNuos+GggLtITokM8xkTAUFmoqecBN9NhQUaA/RIZlhJmMqKNBU9ISb6LOhoEB7iA7JDDMZU0GBpqIn3ESfDQUF2kN0SGaYyZgKCjQVPeEm+mwoKNAeokMyw0zGVFCgqegJN9FnQ0GB9hAdkhlmMqaCAk1FT7iJPhsKCrSH6JDMMJMxFRRoKnrCTfTZUFCgPUSHZIaZjKmgQFPRE26iz4aCAu0hOiQzzGRMBQWaip5wE302FBRoD9EhmWEmYyoo0FT0hJvos6GgQHuIDskMMxlTsWmBfvPW9xp/f/bzOwcHP/h0u4OJnnATPeFmuyyHQYH2EB2SGWYypmLTAv3koFGg37x1UPr232x1MNETbqIn3GwV5UAo0B6iQzLDTMZUbFagzz45aBboJwevfHr69XsHr/xhm4OJnnATPeFmmySHQoH2EB2SGWYypmKjAv39jw+aBfrVneq25zdvvfyX2xxM9ISb6Ak3ZKUjK902G4q2TQr0i4ODH/6uUaBf2N++OPjRNgcTPeEmesINWenISrfNhqJtowL97r86fdgo0E8OflL9b/OtvUVPuImecENWOrLSbbOhaNv0QaRGVT57z+66f3Wn/UPQL/uInnATPeGGrHRkpeu1khvWw3RQoKnoCTdkpSMrXa+V3LAepsO7QLd6IlP0hJvoCTdkpSMr3RYLisSub4H2Ej3hJnrCDVnpyEq3xYIiQYGmoifckJWOrHRbLCgSLgXKo/C7QFY6stJts6Fo8ynQxfM/eR6oI7LSkZVumw1Fm0+B8kqkHSArHVnpttlQtPkU6LP3Dr7La+GdkZWOrHTbbCjati3Qr+5UNzq/5mpM7shKR1a6rVYULU4Fevr1z2f9+YOtbn8y6E1kpSMr3XY7iiauSJ+KnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376CAs2JnnBDVjqy0vntKyjQnOgJN2SlIyud376iq0BP3vn+g+VfHt38Ow+y7+UuesJN9IQbstKRlc5vX9FVoE9f/dav8n/ZqegJN9ETbshKR1Y6v32FVKCP9inQCGSlIyud374iU6BPXy0Sl7kLH4CsdGSl89tX5G6BHqcFenugg4mecBM94YasdGSl89tX5Ar05K9u3bq5v3ft1sLrHw11MNETbqIn3JCVjqx0fvsK6Wegw4mecBM94YasdGSl89tXSE9jGk70hJvoCTdkpSMrnd++gifS50RPuCErHVnp/PYV6wv0twufD3Qw0RNuoifckJWOrHR++4ruAn38du1ReJ4HGoGsdGSl89tXdBZo89mgFGgEstKRlc5vX9FZoEdFcen1Xyz8kifSByArHVnp/PYVnY/C3y2uDHwgpegJN9ETbshKR1Y6v31F9/NA93428IGUoifcRE+4ISsdWen89hU8kT4nesINWenISue3r+i+C88t0HBkpSMrnd++Ys2DSNeHPY5K9ISb6Ak3ZKUjK53fvqKzQGc3Qd8c9kBK0RNuoifckJWOrHR++4ru18LfLIrVBZmGemF89ISb6Ak3ZKUjK53fvqL7QaSCJ9JHIysdWen89hUUaE70hBuy0pGVzm9fwdWY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Kx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRoTvSEG7LSkZXOb19BgeZET7ghKx1Z6fz2FRRo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)

> **Interpreting scenario differences**
>
> These scenario differences are changes in long-run stochastic
> behavior, not forecasts for any single year. A higher `mean_total` or
> `mean_p_hur` means the simulated climate state is more active on
> average, but quiet years still occur.

## 4 Daily Series Comparison

Annual counts are useful for climatology, but many planning questions
depend on daily disruption. To illustrate that step, we generate daily
hazard-impact series for Saba under the baseline and the warmer
scenario.

``` r
daily_baseline <- generate_daily_hazard_impact_spatial(
  out = out_baseline,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = "Baseline",
  seed = 123
)$Saba

daily_future <- generate_daily_hazard_impact_spatial(
  out = out_future,
  location = "Saba",
  sim_years = 1:200,
  year0 = 2025,
  gust_factor = 1.25,
  damage = list(method = "intensity"),
  pulse_shape = "cosine",
  scenario = "Future (+1.5 C)",
  seed = 123
)$Saba
```

The combined table lets us compare how often Saba experiences
tropical-storm and hurricane days, along with the average annual
damage-rate signal over the 200 simulated years.

``` r
daily_summary <- bind_rows(daily_baseline, daily_future) |>
  group_by(scenario) |>
  summarise(
    tc_days = sum(wind_kt > 0, na.rm = TRUE),
    ts_days = sum(wind_kt >= 34, na.rm = TRUE),
    hur_days = sum(wind_kt >= 64, na.rm = TRUE),
    mean_annual_damage = sum(damage_rate, na.rm = TRUE) / 200,
    .groups = "drop"
  )

daily_summary
```

    #> # A tibble: 2 × 5
    #>   scenario        tc_days ts_days hur_days mean_annual_damage
    #>   <chr>             <int>   <int>    <int>              <dbl>
    #> 1 Baseline            840     408       81            0.00314
    #> 2 Future (+1.5 C)     863     436       88            0.00291

``` r
daily_plot_data <- data.frame(
  scenario = rep(daily_summary$scenario, times = 3),
  metric = rep(
    c("TC days/year", "TS days/year", "HUR days/year"),
    each = nrow(daily_summary)
  ),
  value = c(
    daily_summary$tc_days / 200,
    daily_summary$ts_days / 200,
    daily_summary$hur_days / 200
  )
)

daily_plot_data |>
  ggplot(aes(x = metric, y = value, fill = scenario)) +
  geom_col(position = "dodge") +
  labs(
    x = NULL,
    y = "Mean annual days",
    fill = "Scenario",
    title = "Daily disruption metrics for Saba by scenario"
  ) +
  theme_minimal(base_size = 11)
```

![](data:image/png;base64,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)

| Column | Description |
|----|----|
| `tc_days` | Number of simulated days with any tropical-cyclone wind signal (`wind_kt > 0`). |
| `ts_days` | Number of days reaching at least 34 kt. |
| `hur_days` | Number of days reaching at least 64 kt. |
| `mean_annual_damage` | Average annual damage-rate total across the 200 simulated years. |

> **Adding trade-wind background across both scenarios**
>
> The daily series above have zero wind on all non-TC days. For
> continuous-wind applications (energy, chronic exposure, port
> disruption), embed a
> [`make_background_wind_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_background_wind_cfg.md)
> object in each hazard config so it flows through automatically to
> [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md):
>
> ``` r
> bg_cfg <- make_background_wind_cfg(
>   weibull_params = list(
>     Saba = tibble::tibble(month = 1:12, shape = 2.0, scale = 13)
>   ),
>   ar1 = 0.4
> )
>
> cfg_baseline_bg <- make_hazard_cfg(
>   data_path = ibtracs_demo,
>   historical_start_year = 1970L,
>   simulation_years = simulation_years,
>   climate = make_climate_cfg(scenario = "stationary"),
>   background_wind = bg_cfg            # <-- attach once to each cfg
> )
>
> cfg_future_bg <- make_hazard_cfg(
>   data_path = ibtracs_demo,
>   historical_start_year = 1970L,
>   simulation_years = simulation_years,
>   climate = make_climate_cfg(delta_sst = 1.5, sensitivity_mode = "fixed"),
>   background_wind = bg_cfg
> )
> ```
>
> After
> [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md),
> call
> [`generate_daily_hazard_impact_spatial()`](https://tanerumit.github.io/ipdcstorm/reference/generate_daily_hazard_impact_spatial.md)
> without a `background_wind` argument — it reads the setting from
> `out$cfg$background_wind` automatically, keeping the same seed and CRN
> pattern so the climate comparison remains clean. See [Tutorial 1
> §Background
> Wind](https://tanerumit.github.io/ipdcstorm/articles/tutorial_1_setup.qmd#sec-background-wind)
> for parameter guidance and a worked example.

> **About Level 3 perturbation**
>
> The package also supports Level 3 storm perturbation for more
> exploratory scenario analysis. We do not use it here, so the daily
> comparison isolates the Level 1 and Level 2 climate effects discussed
> above.

## 5 Key Parameters

For most climate sensitivity runs, three settings matter most:

| Parameter | What it controls |
|----|----|
| `delta_sst` | The imposed MDR SST shift in °C relative to the 1991-2020 baseline. In this tutorial, `delta_sst = 1.5` defines the future climate state directly. |
| `scenario` vs `delta_sst` | Use a named `scenario` when you want a built-in helper tied to a source such as IPCC AR6. Use `delta_sst` when you want a transparent what-if value without committing to one named pathway. |
| `sensitivity_mode` | `"fixed"` keeps the calibrated climate sensitivities unchanged as warming increases. `"linear_shifted"` lets those sensitivities themselves vary with `delta_sst`, which is more exploratory. |

This baseline-versus-scenario pattern is a good starting point for
stress testing. Once you are comfortable with the annual and daily
outputs, you can expand the same workflow to additional islands or
alternative warming assumptions.
