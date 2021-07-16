"""
Find a reasonable fear function
parameter :
- I
- essentials
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import date

"""

The columns are :
['Date_YMD', 'Daily Confirmed', 'Total Confirmed', 'Daily Recovered', 'Total Recovered',
'Daily Deceased','Total Deceased']

2020-04-15 - 2021-02-09 : First wave
2021-02-10 - Now: Second Wave

2020-02-15 is row index 16
2020-04-15 is row index 76
2020-08-01 is row index 184
2021-02-10 is row index 377
"""
covid_data = pd.read_csv(
    "https://api.covid19india.org/csv/latest/case_time_series.csv"
).drop(["Date"], axis=1)

# convert date from string to python date object
dates = covid_data.values[:, 0]
date_obj = list(map(date.fromisoformat, dates))
covid_data["Date_YMD"] = date_obj

R_1 = covid_data["Total Recovered"] + covid_data["Total Deceased"]
I_1 = covid_data["Total Confirmed"] - R_1


"""
The columns are :
['date', 'retail_and_recreation_percent_change', 'grocery_and_pharmacy_percent_change',
'parks_percent_change','transit_stations_percent_change', 'workplaces_percent_change',
'residential_percent_change', 'sum_except_groceries_and_pharmacy']

data starts from 2020-02-15
2020-08-01 is index 168
"""
mobility_data = pd.read_csv("India_COVID_Mobility.csv")

# convert date from string to python date object
dates = mobility_data.values[:, 0]
processed_dates = [list(map(int, day.split("/"))) for day in dates]
date_obj = [date(yr, mnth, day) for (mnth, day, yr) in processed_dates]
mobility_data["date"] = date_obj


# plot Active Cases & 'Fear' side by side vs time since 2020-08-01

I_range = I_1[184:]
fear = mobility_data["sum_except_groceries_and_pharmacy"][168:]
essensial = -mobility_data["grocery_and_pharmacy_percent_change"][168:]

fig, ax1 = plt.subplots()

color = "tab:red"
ax1.set_xlabel("time")
ax1.set_ylabel("Active Cases", color=color)
ax1.plot_date(
    covid_data["Date_YMD"][184:], I_range, "-", label="Active Cases", color=color,
)
ax1.tick_params(axis="y", labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = "tab:blue"
ax2.set_ylabel(
    "sum except groceries and pharmacy", color=color
)  # we already handled the x-label with ax1
ax2.plot_date(
    mobility_data["date"][168:],
    fear,
    "-",
    label="sum_except_groceries_and_pharmacy",
    color=color,
)
ax2.tick_params(axis="y", labelcolor=color)

ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = "tab:green"
ax3.set_ylabel(
    "grocery and pharmacy percent change", color=color
)  # we already handled the x-label with ax1
ax3.plot_date(
    mobility_data["date"][168:],
    essensial,
    "-",
    label="grocery and pharmacy percent change",
    color=color,
)
ax3.tick_params(axis="y", labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.legend()


# plot Active cases vs 'fear' between 2021-02-10 2021-07-07

# I_range = I_1[377:525]
# maxI = np.max(I_range)

# mobility = mobility_data["sum_except_groceries_and_pharmacy"][361:]
# maxmobility = np.max(mobility)
# plt.plot(mobility / maxmobility, I_range / maxI, ".")


plt.show()
