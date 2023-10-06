from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import os
import prince
import matplotlib.pyplot as plt

pa = "/home/kshedden/data/Teaching/bhht"

pdf = PdfPages("bhht_py_mca.pdf")

df = pd.read_csv(os.path.join(pa, "cross-verified-database.csv.gz"), encoding="latin-1")

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

cols = {"occ": "orange", "sex": "purple", "reg": "lime", "bcen": "navy"}

# Plot the category scores.  There are so many object scores that a scatterplot
# of object scores would be difficult to interpret.
for (x, y) in [(0, 1), (0, 2), (1, 2)]:
    cc = mca.column_coordinates(df)

    xmin, xmax = cc.iloc[:, x].min(), cc.iloc[:, x].max()
    d = xmax - xmin
    xmin -= 0.1*d
    xmax += 0.1*d
    ymin, ymax = cc.iloc[:, y].min(), cc.iloc[:, y].max()
    d = ymax - ymin
    ymin -= 0.1*d
    ymax += 0.1*d

    plt.clf()
    plt.grid(True)
    for k in cols.keys():
        cx = cc[cc.index.str.startswith(k)]
        if k == "bcen":
            plt.plot(cx.iloc[:, x], cx.iloc[:, y], "-", color=cols[k])
        for i in range(cx.shape[0]):
            plt.text(cx.iloc[i, x], cx.iloc[i, y], cx.index[i], color=cols[k],
                     ha="center", va="center")
    plt.xlabel("Component %d" % (x + 1))
    plt.ylabel("Component %d" % (y + 1))
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    pdf.savefig()

pdf.close()
