from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import prince
import matplotlib.pyplot as plt

pdf = PdfPages("bhht_py_mca.pdf")

df = pd.read_csv("cross-verified-database.csv.gz", encoding="latin-1")

df = df.rename({"level1_main_occ": "occ", "gender": "sex", "un_region": "reg"}, axis=1)
df = df.loc[:, ["birth", "occ", "sex", "reg"]]
df = df.dropna()
df = df.loc[df.birth >= 1500, :]

# Remove very infrequent or difficult to interpret categories
df = df.loc[df.occ != "Other", :]
df = df.loc[df.occ != "Missing", :]
df = df.loc[df.sex != "Other", :]

# Create a century of birth variable, make it a string so that
# it is interpreted as a nominal variable.
df.loc[:, "bcen"] = df.birth.round(-2)
df = df.drop("birth", axis=1)
df.loc[:, "bcen"] = ["%d" % x for x in df.bcen]

mca = prince.MCA(n_components=3)
mca = mca.fit(df)
mca.transform(df)

# Plot the category scores.  There are so many object scores that a naive scatterplot
# of object scores would be difficult to interpret.
for (x, y) in [(0, 1), (0, 2), (1, 2)]:
    ax = mca.plot_coordinates(X=df, x_component=x, y_component=y, show_row_points=False, show_column_labels=True)
    pdf.savefig(ax.get_figure())

pdf.close()
