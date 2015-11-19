//
//  main.cpp
//  RevisedSimplex
//
//  Created by Xuan Baihua on 10/18/15.
//  Copyright (c) 2015 Xuan Baihua. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

// zero tolerances
static const double epsilon1 = 0.00001;
static const double epsilon2 = 0.00000001;

// number of constraints
static size_t m;

// number of variables
static size_t n;

// a structure that represents a basic variable
// containing its label (ranges from 0 to m + n - 1) and its current value
struct variable {
    size_t label;
    double value;
};

// a structure that represents an eta matrix
// containing the column that's different from an identity matrix (ranges from 1 to m - 1)
// and the values in that column
struct eta {
    size_t col;
    vector<double> values;
};

// print out the A matrix
void printMatrix(double *matrix, size_t r, size_t c) {
    for (size_t row = 0; row < r; ++row) {
        for (size_t col = 0; col < c; ++col) {
            printf("%10.3f ", matrix[row * (m + n) + col]);
        }
        printf("\n");
    }
}

// print out the information of the LP, which includes
// number of variables, constraints, coefficients in the objective function,
// initial values of the basic variables and the A matrix
void printLPInfo(double objFuncCoeff[], variable b[], double *matrix) {
    printf("m = %lu ", m);
    printf("n = %lu \n",n);
    
    printf("c = ");
    
    for (size_t i = 0; i < m + n; ++i) {
        printf("%10.3f ", objFuncCoeff[i]);
    }
    
    printf("\nb = ");

    for (size_t i = 0; i < m; ++i) {
        printf("%10.3f ", b[i].value);
    }

    printf("\nA = \n");
    
    printMatrix(matrix, m, n + m);
};

// print out basic and non-basic variables for each iteration
void printVariables(size_t nonbasic[], variable b[]) {
    printf("N = { ");
    
    for (size_t i = 0; i < n; ++i) {
        printf("x%lu ", nonbasic[i] + 1);
    }
    
    printf("} B = { ");
    
    // basic variables
    for (size_t i = 0; i < m; ++i) {
        printf("x%lu ", b[i].label + 1);
    }
    
    printf("}\n");
};

// print out the values of basic variables for each iteration
void printBbar(variable b[]) {
    printf("bbar = ");
    
    for (size_t i = 0; i < m; ++i) {
        printf("%10.3f ", b[i].value);
    }
    
    printf("\n");
}

// print out the final values of decision variables and slack variables
// when the optimal solution is reached
void printFinalVariables(variable b[], size_t nonbasic[]) {
    double varValues[m + n];
    
    for (size_t row = 0; row < m; ++row) {
        varValues[b[row].label] = b[row].value;
    }
    
    for (size_t col = 0; col < n; ++col) {
        varValues[nonbasic[col]] = 0.0;
    }
    
    // decision variables
    printf("Decision variables: ");
    
    for (size_t i = 0; i < n; ++i) {
        printf("x%lu = %5.3f ", i + 1, varValues[i]);
    }
    
    // slack variables
    printf("\nSlack variables: ");
    
    for (size_t i = n; i < m + n; ++i) {
        printf("x%lu = %5.3f ", i + 1, varValues[i]);
    }
    
    printf("\n");
}

// function to print out the family of solutions when the LP is found to be unbounded
void printFamilyOfSolutions(variable b[], size_t nonbasic[], vector<double> d, double largestCoeff, size_t enteringLabel, double z) {
    // in order to print out variable information in their natural order
    // we need to remember which row a basic variable locates in the basis
    // as to correspond it to the correct row in the entering variable's column in the dictionary
    struct varInfo {
        size_t label;
        double value;
        size_t row_in_basis;
        bool isBasic;
    };
    
    varInfo info[m + n];
    
    for (size_t row = 0; row < m; ++row) {
        variable v = b[row];
        info[v.label] = {v.label, v.value, row, true};
    }
    
    for (size_t col = 0; col < n; ++col) {
        info[nonbasic[col]] = {nonbasic[col], 0.0, 0, false};
    }
    
    // decision variables
    printf("Decision variables: ");
    
    for (size_t i = 0; i < n; ++i) {
        printf("x%lu = %5.3f ", i + 1, info[i].value);
        if (info[i].isBasic) {
            printf("+ %5.3fx%lu ", d[info[i].row_in_basis] * (-1.0), enteringLabel + 1);
        }
    }
    
    // slack variables
    printf("\nSlack variables: ");
    
    for (size_t i = n; i < m + n; ++i) {
        printf("x%lu = %5.3f ", i + 1, info[i].value);
        if (info[i].isBasic) {
            printf("+ %5.3fx%lu ", d[info[i].row_in_basis] * (-1.0), enteringLabel + 1);
        }
    }
    
    // objective function
    printf("\nZ = %5.3f + %5.3fx%lu, ", z, largestCoeff, enteringLabel + 1);
    
    printf("with x%lu >= 0\n", enteringLabel + 1);
}


// a function for sorting variables in a vector
// into descending order
bool mComparator(variable v1, variable v2) {
    return v1.value > v2.value;
}

int main(int argc, const char * argv[]) {
    
    // Read in the number of constraints and
    // the number of variables
    cin >> m >> n;
    
    // Coefficients in the objective function
    double objFuncCoeff[n + m];
    
    // Read in the coefficients in the objective function
    for (size_t col = 0; col < n; ++col) {
        cin >> objFuncCoeff[col];
    }
    
    // Populate the rest of C to be 0
    for (size_t col = n; col < m + n; ++col) {
        objFuncCoeff[col] = 0.0;
    }
    
    // Matrix A represented as a 1-D array
    // with rows stacked next to each other
    double A[m * (n + m)];
    
    // Column of labels and values of the basic variables in the basic feasible solution
    variable b[m];
    
    // Column of variable (only with labels, from 0 to m + n - 1) of the nonbasic variables
    size_t nonbasic[n];
    
    // Initialize columns for decision variables in matrix A
    // Initialize b, the array of basic variables
    for (size_t row = 0; row < m; ++row) {
        for (size_t col = 0; col <= n; ++col) {
            if (col == n) {
                double bRow;
                cin >> bRow;
                variable bVar = {n + row, bRow};
                b[row] = bVar;
            } else {
                cin >> A[row * (m + n) + col];
            }
        }
    }
    
    // Initialize columns for slack variables in matrix A
    for (size_t row = 0; row < m; ++row) {
        size_t base = (m + n) * row + n;
        for (size_t col = 0; col < m; ++col) {
            if (col != row) {
                A[base + col] = 0.0;
            } else {
                A[base + col] = 1.0;
            }
        }
    }
    
    // Initialize the nonbasic variables' labels to be 0 through (n - 1)
    for (size_t i = 0; i < n; ++i) {
        nonbasic[i] = i;
    }
    
    // Print out initial input values
    printLPInfo(objFuncCoeff, b, A);
    
    printf("\n\n");
    
    // Print out initial nonbasic and basic variables
    printVariables(nonbasic, b);
    
    // Print out initial values of basic variables
    printBbar(b);
    
    printf("\n");
    
    // Check initial feasibility
    for (size_t row = 0; row < m; ++row) {
        if (b[row].value < 0.0) {
            printf("The given linear program is infeasible, exiting the program.\n");
            return 0;
        }
    }
    
    // Initial basic solution is feasible, now proceed with the Simplex Method
    
    // A counter to remember the current iteration number
    size_t counter = 1;
    
    // An array of eta matrices representing previous pivots
    vector<eta> pivots{};
    
    // Initial value of objective function
    double z = 0.0;
    
    // Revised Simplex Method
    while (true) {
        printf("Iteration%lu\n------------\n", counter);
        
        // compute y using eta matrices (yB = Cb)
        vector<double> y(m);
        
        // initialize y to be Cb
        for (size_t row = 0; row < m; ++row) {
            variable v = b[row];
            y[row] = objFuncCoeff[v.label];
        }
        
        // solving y in yB = Cb
        for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter) {
            eta pivot = *rIter;
            size_t colToChange = pivot.col;
            double yOriginal = y[colToChange];
            
            for (size_t row = 0; row < pivot.values.size(); ++row) {
                if (row != colToChange) {
                    yOriginal -= pivot.values[row] * y[row];
                }
            }
            
            double yNew = yOriginal / pivot.values[colToChange];
            y[colToChange] = yNew;
        }
        
        // print out solved y
        printf("y = ");
        
        for (auto iter = y.cbegin(); iter != y.cend(); ++iter) {
            printf("%10.3f ", *iter);
        }
        
        printf("\n");
        
        // choose an entering column
        // using the condition Cn > ya, where "a" is a column of An
        
        // a vector to keep track of the variables
        // whose coefficients in the objective function in this iteration are positive
        vector<variable> cnbars;
        
        size_t enteringLabel = nonbasic[0];
        double largestCoeff = -1.0;
        
        // print cnbar
        printf("cnbar: ");
        
        for (size_t i = 0; i < n; ++i) {
            size_t varLabel = nonbasic[i];
            double cni = objFuncCoeff[varLabel];
            double yai = 0.0;
            
            for (size_t yIndex = 0; yIndex < m; ++yIndex) {
                yai += y[yIndex] * A[yIndex * (m + n) + varLabel];
            }
            
            double cnbar = cni - yai;
            
            printf("x%lu %5.3f ", varLabel + 1, cnbar);
            
            if (cnbar > epsilon1) {
                variable v = {varLabel, cnbar};
                
                cnbars.push_back(v);
                
                if (cnbar > largestCoeff) {
                    largestCoeff = cnbar;
                    enteringLabel = varLabel;
                }
            }
        }
        
        // sort the variables into descending order
        // based on their coefficients in the objective function
        sort(cnbars.begin(), cnbars.end(), mComparator);
        
        printf("\n");
        
        // If the vector cnbars is empty, then there are no candidates for the entering variable
        
        if (cnbars.size() == 0) {
            printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
            printFinalVariables(b, nonbasic);
            return 0;
        } else {
            printf("Entering variable is x%lu \n", enteringLabel + 1);
        }
        
        size_t enteringVariable_index = 0;
        
        // compute the column d in Anbar
        // for the entering variable
        // using eta matrices (Bd = a)
        vector<double> d(m);
        
        size_t leavingLabel;
        size_t leavingRow;
        double smallest_t;
        
        while (true) {
            
            leavingLabel = -1;
            leavingRow = -1;
            smallest_t = -1;
            
            if (enteringVariable_index > 0) {
                printf("\n\nRechoosing entering variable since the diagonal element in the eta column is close to zero.\n");
            }
            
            if (enteringVariable_index < cnbars.size()) {
                enteringLabel = cnbars[enteringVariable_index].label;
                
                if (enteringVariable_index > 0) {
                    printf("Entering variable is x%lu \n", enteringLabel + 1);
                }
            } else {
                printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
                printFinalVariables(b, nonbasic);
                return 0;
            }
            
            // initialize d to be the entering column a
            for (size_t row = 0; row < m; ++row) {
                d[row] = A[row * (m + n) + enteringLabel];
            }
            
            // Go through eta matrices from pivot 1 to pivot k
            for (auto iter = pivots.cbegin(); iter != pivots.cend(); ++iter) {
                eta pivot = *iter;
                size_t rowToChange = pivot.col;
                double dOriginal = d[rowToChange];
                
                d[rowToChange] = dOriginal / pivot.values[rowToChange];
                
                for (size_t row = 0; row < d.size(); ++row) {
                    if (row != rowToChange) {
                        d[row] = d[row] - pivot.values[row] * d[rowToChange];
                    }
                }
            }
            
            // print out d (abarj)
            printf("d = ");
            
            for (auto iter = d.cbegin(); iter != d.cend(); ++iter) {
                printf("%5.3f ", *iter);
            }
            
            printf("\n");
            
            // compute t for each b[i].value / d[i]
            // where d[i] > 0
            // choose the corresponding i for the smallest ratio
            // as the leaving variable
            
            // initialize smallest_t to be the first ratio where
            // the coefficient of the entering variable in that row is negative
            for (size_t row = 0; row < d.size(); ++row) {
                if (d[row] > 0.0) {
                    leavingLabel = b[row].label;
                    leavingRow = row;
                    smallest_t = b[row].value / d[row];
                }
            }
            
            // if no ratio is computed, then the LP is unbounded
            if (leavingLabel == -1) {
                printf("\nThe given LP is unbounded. The family of solutions is:\n");
                printFamilyOfSolutions(b, nonbasic, d, largestCoeff, enteringLabel, z);
                return 0;
            }
            
            // there is at least one ratio computed, print out the ratio(s)
            // and choose the row corresponding to the smallest ratio to leave
            printf("ratio: ");
            
            for (size_t row = 0; row < d.size(); ++row) {
                if (d[row] < 0.0) {
                    continue;
                }
                
                double t_row = b[row].value / d[row];
                
                if (t_row >= 0.0) {
                    printf("x%lu %5.3f ", b[row].label + 1, t_row);
                }
                
                if (t_row < smallest_t) {
                    leavingLabel = b[row].label;
                    leavingRow = row;
                    smallest_t = t_row;
                }
            }

            // check the diagonal element in the eta column
            // to see if the current choice of entering variable has to be rejected
            if (d[leavingRow] > epsilon2) {
                printf("\nLeaving variable is x%lu\n", leavingLabel + 1);
                break;
            } else {
                enteringVariable_index++;
                continue;
            }
        }
        
        // At this point we have a pair of entering and leaving variables
        // so that the entering variable is positive and the diagonal entry in the eta column
        // of the eta matrix is fairly far from zero
        
        // set the value of the entering varaible at t
        // modify b (change leaving variable to entering variable, change values of other basic vars)
        variable enteringVar = {enteringLabel, smallest_t};
        b[leavingRow] = enteringVar;
        
        for (size_t row = 0; row < sizeof(b) / sizeof(b[0]); ++row) {
            if (row != leavingRow) {
                b[row].value -= d[row] * smallest_t;
            }
        }
        
        // push a new eta matrix onto the vector
        eta pivot = {leavingRow, d};
        pivots.push_back(pivot);
        
        // print out the eta matrix representing the pivot at this iteration
        printf("E%lu = column %lu: ", counter, leavingRow);
        
        for (auto iter = d.cbegin(); iter != d.cend(); ++iter) {
            printf("%5.3f ", *iter);
        }
        
        printf("\n");
        
        nonbasic[enteringLabel] = leavingLabel;
        
        // print out nonbasic and basic variable set after the pivot
        printVariables(nonbasic, b);
        
        // print out the new values of the basic variables
        printBbar(b);
        
        // print out the coefficient of the entering variable and the amount the entering variable has been increased
        printf("\nCoefficient of entering variable: %5.3f\nAmount increased for the entering variable is: %5.3f\n", largestCoeff, smallest_t);
        
        // increase the value of the objective function
        double increasedValue =largestCoeff * smallest_t;
        
        // print out the update to the objective function value
        printf("Increased value: %5.3f\n", increasedValue);
        
        double originalZ = z;
        
        z += increasedValue;
        
        printf("Value of the objective function changed from %5.3f to %5.3f\n\n\n", originalZ, z);
        
        counter++;
    }
    
    return 0;
}