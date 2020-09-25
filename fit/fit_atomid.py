#MIT License
#
#Copyright (c) 2020 TJ O'Donnell
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
import pandas as pd
from sklearn.linear_model import LinearRegression, BayesianRidge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
from sklearn.metrics import r2_score
import sqlite3
import pickle
import sys
import os
import json
import time
import matplotlib.pyplot as plt

def find_correlation(data, threshold=0.9, remove_negative=False):
    import numpy as np
    """
    https://gist.github.com/Swarchal/e29a3a1113403710b6850590641f046c

    modified tjo 2020-21-July
    . to return a list of correlated items instead of just flat list of items to remove
    . to check for correlation values >= instead of >

    Given a numeric pd.DataFrame, this will find highly correlated features,
    and return a list of features to remove.
    Parameters
    -----------
    data : pandas DataFrame
        DataFrame
    threshold : float
        correlation threshold, will remove one of pairs of features with a
        correlation greater than this value.
    remove_negative: Boolean
        If true then features which are highly negatively correlated will
        also be returned for removal.
    Returns
    --------
    select_flat : list
        listof column names to be removed
    """
    corr_mat = data.corr()
    if remove_negative:
        corr_mat = np.abs(corr_mat)
    corr_mat.loc[:, :] = np.tril(corr_mat, k=-1)
    already_in = set()
    result = []
    for col in corr_mat:
        perfect_corr = corr_mat[col][corr_mat[col] >= threshold].index.tolist()
        #corr_values = corr_mat[col][corr_mat[col] >= threshold].values.tolist()
        if perfect_corr and col not in already_in:
            already_in.update(set(perfect_corr))
            perfect_corr.append(col)
            result.append(perfect_corr)
            #result.append(corr_values)
    return result
#    select_nested = [f[1:] for f in result]
#    return select_nested
#    select_flat = [i for j in select_nested for i in j]
#    return select_flat

def parentage(con, atomid):
    cur = con.cursor()
    result = []
    sql = """With Recursive parent_of(atomid, parentid, iteration) As (
              Select atomid, parentid, iteration From parents Where atomid = ?
             Union
              Select parents.atomid, parents.parentid, parents.iteration From parents, parent_of Where parents.atomid = parent_of.parentid
             )
             Select atomid, parentid, iteration, smarts From parent_of Left Join smarts Using (atomid) Order By iteration Desc"""
    cur.execute(sql, [atomid])
    for row in cur:
        result.append({"atomid":row[0], "parentid":row[1], "iteration":row[2], "smarts":row[3]})
    return result

def getCounts(cur, nmols, max_level, flimit):
    # get high-frequency atomid to use as features
    cur.execute("Create Temporary Table If Not Exists top As Select * From atom_types Where atomid Is Not Null And iteration <= %d And frequency > 1 Order By (iteration=0) Desc, frequency Desc Limit %d" % (max_level, flimit))
    #cur.execute("Create Temporary Table If Not Exists top As Select * From atom_types Where iteration Between 1 and %d Order By frequency Desc Limit %d" % (max_level, flimit))
    cur.execute("Select atomid, iteration, frequency From top Order By atomid")
    atomid = {"atomid": [], "iteration": [], "frequency": []}
    for row in cur:
        atomid["atomid"].append(row[0])
        atomid["iteration"].append(row[1])
        atomid["frequency"].append(row[2])
    natomid = len(atomid["atomid"])
    #print (natomid)

    # get counts of occurrences of each atomid for each molecule
    cur.execute("Create Temporary Table If Not Exists counts As Select molid, atomid, Count(atom_index) pcount From atoms Group By molid, atomid")
    cur.execute("With matrix As (Select atomid, molid From top Join fitprop) Select Coalesce(pcount,0) pcount From matrix Left Join counts Using (molid,atomid) Order By molid, atomid")
    counts = []
    for i in range(0, nmols):
        counts.append([row[0] for row in cur.fetchmany(natomid)])
    #print (len(counts), len(counts[0]))
    #print (counts)
    return (atomid, counts)

def printElapsedTime(start_time):
    # print elapsed time since last call
    end_time = time.time()
    if start_time: print ("%d seconds." % (end_time - start_time))
    return end_time
    
def order_atomid(curr, c):
    ord = []
    cur.execute("Select atomid From atom_types Where atomid in (%s) Order By iteration" % ','.join(['?']*len(c)), c)
    for row in cur:
        ord.append(row[0])
    return ord

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Create a model that fits molecular properties to atomid counts",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("database", help="input structures+atomid+properties sqlite3 (db3) file")
    parser.add_argument("property", help="name of property (sdf tag) to fit")
    parser.add_argument("model", help="output model sqlite3 (db3) file")
    parser.add_argument("ndescriptors", type=int, help="number of atomid descriptors to consider")
    parser.add_argument("level", type=int, help="maximum atomid level to consider")
    parser.add_argument("-p", "--pickle", help="output file for python pickle of model", default=None)
    parser.add_argument("-l", "--list", help="list properties in input db3 file, and exit", action="store_true")
    parser.add_argument("-f", "--fit", help="method of fitting", choices=["mlr", "bayes", "rf", "xrt"], default="mlr")
    parser.add_argument("-c", "--correlated", help="keep or remove correlated atomid before fitting", choices=["remove", "reChcek", "keep"], default="remove")
    parser.add_argument("-n", "--modulo_N", help="include only every nth molecule; useful for shorter tests", type=int, default=1)
    parser.add_argument("--plot", help="output file for plot of input vs predicted values", default=None)
    parser.add_argument("-v", "--verbosity", type=int, help="increase output verbosity", default=1)
    return parser

def main():
    t = printElapsedTime(0)
    parser = parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    parsed = parser.parse_args()
    
    mol_db = parsed.database
    property_name = parsed.property
    model_db = parsed.model
    flimit = parsed.ndescriptors
    max_level = parsed.level
    method = parsed.fit
    verbosity = parsed.verbosity
    if parsed.correlated == "keep":
        remove_correlated = 0
    elif parsed.correlated == "remove":
        remove_correlated = 1
    elif parsed.correlated == "reCheck":
        remove_correlated = 2
    plotfile = parsed.plot
    mod_n = parsed.modulo_N
    list_prop = parsed.list
    pickle_file = parsed.pickle

    if not os.path.isfile(mol_db):
        print("cannot open file %s" % mol_db)
        sys.exit()
    con = sqlite3.connect(mol_db)
    cur = con.cursor()
    
    # list properties and exit
    if list_prop:
        cur.execute("Select propname From property_names")
        for row in cur:
            print (row[0])
        sys.exit()
    
    # don't overwrite output file.
    if os.path.isfile(model_db):
        print ("cannot overwrite existing file %s" % model_db)
        sys.exit()
     # don't overwrite output file
    if pickle_file and os.path.isfile(pickle_file):
        print ("cannot overwrite existing file %s" % pickle_file)
        sys.exit()
    
    # model output database
    mcon = sqlite3.connect(model_db)
    mcur = mcon.cursor()
    mcur.execute("Create Table If Not Exists parameters (start_date Text, input_file Text, property_name Text, model_file Text, n_features Integer, max_level Integer, method Text)")
    mcur.execute("Create Table If Not Exists results (end_date Text, score Numeric, singular Integer)")
    mcur.execute("Insert Into parameters Values (datetime('now'),?,?,?,?,?,?)", (mol_db, property_name, model_db, flimit, max_level, method))

    # get list of property values
    sql = """Create Temporary Table fitprop As
        With vtemp As (Select molid, Cast(propvalue As Float) As value From property_values
        Join property_names Using (propid) Where propname='%s' And propvalue Not Like '<%%' And propvalue Not Like '>%%')
        Select molid, vtemp.value From molecule Left Join vtemp Using (molid) Where vtemp.value Is Not Null And molid %% %d = 0"""
    cur.execute(sql % (property_name, mod_n))
    
    property_values = []
    cur.execute("Select value From fitprop Order By molid")
    for row in cur:
        property_values.append(row[0])
    nmols = len(property_values)
    if verbosity > 0: print ("%d mol properties" % nmols)
    if verbosity > 1: t = printElapsedTime(t)
    if nmols == 0:
        print ("no values for property '%s'" % property_name)
        sys.exit()
    
    (atomid, counts) = getCounts(cur, nmols, max_level, flimit)
    natomid = len(atomid["atomid"])
    if verbosity > 0: print ("%d atomid counted" % natomid)
    if verbosity > 1: t = printElapsedTime(t)
    
    df = pd.DataFrame(counts, columns = atomid["atomid"])
    #print (df.shape)
    correlated = find_correlation(df, threshold = 0.98, remove_negative = True)
    if verbosity > 0: print ("%d correlated atomid sets" % len(correlated))
    if verbosity > 1: t = printElapsedTime(t)
    
    if len(correlated) > 0 and remove_correlated > 0:
        # remove correlated atomid ...
        nremoved = 0
        for c in correlated:
            #print ([[p["atomid"] for p in parentage(con,cc)] for cc in c])
            #print( [cc for cc in c] )
    #         print (".".join(c))
    #         print ()
            #print (order_atomid(cur, c))
            for cc in c[1:]:
                cur.execute("Delete From top Where atomid=?", [cc])
                nremoved += 1
                #print (json.dumps(parentage(con, cc), indent=3))
        if verbosity > 0: print ("%d correlated atomid ignored" % nremoved)
        
        # ... and re-do the count with slimmer set of atomid ...
        (atomid, counts) = getCounts(cur, nmols, max_level, flimit)
        natomid = len(atomid["atomid"])
        if verbosity > 0: print ("%d atomid counted" % natomid)
        if verbosity > 1: t = printElapsedTime(t)
        # ... and look for correlation in reduced set.  should be none
        df = pd.DataFrame(counts, columns = atomid["atomid"])
        if remove_correlated > 1:
            #print (df.shape)
            after_correlated = find_correlation(df, threshold = 0.98, remove_negative = True)
            if verbosity > 0: print ("%d correlated atomid sets" % len(after_correlated))
            if verbosity > 1: t = printElapsedTime(t)
    
    if method == "bayes":
        model = BayesianRidge()
        model.fit(counts, property_values)
        nsingular = len(model.coef_) # ??
        intercept = model.intercept_
        coefficients = model.coef_
    #     score = model.score(counts, property_values)
    #     coefficients = model.coef_
    #     singular = model.coef_
    #     intercept = model.intercept_
    #     print ("Score: %.3f; Intercept: %.3f; %d features" % (score, intercept, len(coefficients)))
    #     property_pred = model.predict(counts)
    elif method == "rf":
        model = RandomForestRegressor(n_estimators=100, min_samples_leaf=1, max_depth = None)
        model.fit(counts, property_values)
        intercept = 0
        coefficients = [None] * natomid
        nsingular = natomid # ??
    elif method == "xrt":
        model = ExtraTreesRegressor(n_estimators=100)
        model.fit(counts, property_values)
        intercept = 0
        coefficients = [None] * natomid
        nsingular = natomid # ??
    else:
        model = LinearRegression().fit(counts, property_values)
        nsingular = len(model.singular_)
        intercept = model.intercept_
        coefficients = model.coef_
    score = model.score(counts, property_values)
    
    if pickle_file:
        pickle.dump(model, open(pickle_file, 'wb'))
        
    #print (len(counts), len(counts[0]))
    property_pred = model.predict(counts)
    r2 = r2_score(property_values, property_pred)
    if verbosity > 0: print ("Score: %.3f; R-squared: %.3f; Intercept: %.3f; %d features; %d singular" % (score, r2, intercept, len(coefficients), nsingular))
    
    # print (len(atomid["atomid"]))
    # print (len(coefficients))
    # for i in range(len(coefficients)):
    #     print (i, atomid["atomid"][i], coefficients[i])
    if verbosity > 1: t = printElapsedTime(t)
    #print (["{0:0.2f}".format(v) for v in property_pred[:20]])
    #print (["{0:0.2f}".format(v) for v in property_values[:20]])
    #print (coefficients)
    if plotfile:
        fig, ax = plt.subplots()
        ax.scatter(property_values, property_pred, marker='.', s=16)
        plt.title("%s %d/%d %s" %(property_name, flimit, max_level, method))
        plt.xlabel("training values")
        plt.ylabel("predicted values")
        fig.savefig(plotfile)
        #plt.show()
    # save coefficients
    mcur.execute("Create Table If Not Exists atomid_coefficients (atomid Text, coefficient Real)")
    mcur.execute("Create Table If Not Exists atomid_info (atomid Text, iteration Integer, frequency Integer)")
    mcur.execute("Insert Into atomid_coefficients (atomid, coefficient) Values (?,?)", ['Intercept', intercept])
    for i in range(0, natomid):
        mcur.execute("Insert Into atomid_coefficients (atomid, coefficient) Values (?,?)", [atomid["atomid"][i], coefficients[i]])
        mcur.execute("Insert Into atomid_info (atomid, iteration, frequency) Values (?,?,?)", [atomid["atomid"][i], atomid["iteration"][i], atomid["frequency"][i]])
    
    # store correlated atomid coefficients
    
    mcur.execute("Create Table correlated (a_atomid Text, b_atomid Text)")
    for c in correlated:
        #mcur.execute("Select coefficient From atomid_coefficients Where atomid In (%s)" % ','.join(['?']*len(c)), c)
        #print ( [ (int(cc), mcur.fetchone()[0]) for cc in c ] )
        a_atomid = c[0]
        for cc in c[1:]:
            mcur.execute("Insert Into correlated (a_atomid, b_atomid) Values (?,?)", (a_atomid, cc))
    
    mcur.execute("Insert Into results Values (datetime('now'),?,?)", (score, nsingular))
    
    mcon.commit()

if __name__ == "__main__":
    import cProfile
    cProfile.run('main()', 'profile.out')
    #main()

