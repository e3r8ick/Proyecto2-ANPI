/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>

#include "Matrix.hpp"
#include "Exception.hpp"

#ifndef ANPI_RESISTORGRID_HPP
#define ANPI_RESISTORGRID_HPP

namespace anpi
{

    ///Pack a pair of indices of the nodes of a resistor
    struct indexPair {
        ///Row of the first node
        std::size_t row1;
        //Column of the first node 
        std::size_t col1;
        ///Row of the second node
        std::size_t row2;
        //Column of the second node 
        std::size_t col2;
    };

class ResistorGrid{
    private:
        ///Matrix of the current equation system
        Matrix<float> A_;
        ///Vector of the current equation system 
        std::vector<float> b_;

        /// Raw map data 
        Matrix<float> rawMap_;
    public:
    
    ///Conastruct
    ResistorGrid(Matrix<float> A, std::vector<float> b, Matrix<float> rawMap){
        A_ = A;
        b_ = b;
        rawMap_ = rawMap;
    } 

    /**
     * Map the index of the terminal nodes
     * */
    std::size_t nodeToIndex(const std::size_t row1,
                            const std::size_t col1,
                            const std::size_t row2,
                            const std::size_t col2)
    {
        if(row1 <= 0 || col1 <= 0 || row2 <= 0 || col2 <= 0){
            throw anpi::Exception("Indices <= 0; no son validos");
        }
        else if (row1 == row2 && col1 == col2)
        {
            throw anpi::Exception("Indices iguales");
        }
        else if (!((((row1+1) == row2 && col1 == col2) || ((row1-1) == row2 && col1 == col2)) || 
                ((row1 == row2 && (col1+1) == col2) || (row1 == row2 && (col1-10) == col2))))
        {
            throw anpi::Exception("Indices no son adyacentes");
        }
        else{
             size_t rowMayor, rowMenor, cooMenor;

            //sacamos el mayor y menor de rows 
            if(row1 >= row2){
                rowMayor = row1;
                rowMenor = row2;
            }else{
                rowMayor = row2;
                rowMenor = row1;
            }

            //sacamos el menor de cols
            if(col1 >= col2){
                cooMenor = col2;
            }else{
                cooMenor = col1;
            }

            size_t n = A_.cols();
            size_t result = 2 * rowMayor * (n-1) + cooMenor + rowMenor + 1;

            if(rowMayor == rowMenor){
                return result;

            }
            return result-(n-1);

        }
    }   

    /// Convert an inde x to the pair of node coordinates
    indexPair indexToNodes (const std::size_t idx){
        anpi::indexPair coordenadas;
        size_t iMayor, iMenor, jMayor, jMenor;
        size_t n = A_.cols();
        size_t mn = 2*n-1;
        size_t mod = idx % mn;
        size_t divi = idx/mn;

        //sacamos j e i mayor y menor
        if (mod == 0){ // last column
            iMayor = divi;
            iMenor = divi - 1;
            jMayor = jMenor = n-1;
            coordenadas = {
                iMenor,
                jMayor,
                jMayor,
                jMenor
            };
        }
        else if(mod > n-1){ //vertical
            mod -= n;
            jMenor = jMayor = mod;
            iMenor = divi;
            iMayor = divi + 1;
            coordenadas = {
                iMenor,
                jMayor,
                iMayor,
                jMenor
            };
        
        }else{ //horizontal
            jMenor = mod -1;
            jMayor = mod;
            iMenor = iMayor = divi;
            coordenadas = {
                iMayor,
                jMenor,
                iMenor,
                jMayor
            };
        }
        return coordenadas;
    }

    //verifica si el nodo es una esquina
    bool isCorner(size_t i, size_t j, size_t m, size_t n){
        if (i == 0 && j == 0){
            return true;
        }else if(i == m-1 && j == n-1){
            return true;
        }else if(i == 0 && j == n-1){
            return true;
        }else if (i == m-1 && j == 0){
            return true;
        }else{
            return false;
        }
    }

    //ecuaciones de nodos
    template <typename T>
    void nodos(size_t m, size_t n, Matrix<T> &A, size_t iin, size_t jin, 
            size_t ifin, size_t jfin, std::vector<T> b){
        bool eliminado = false;
        size_t index = 0;
        size_t m1, m2, m3, m4;
        for(int i=0; i < m; ++i){
            for(int j=0; j<n; ++j){
                if (index >= (2*m*n - m - n)){
                    std::cout << "here" << std::endl;
                }
                
                m1=nodeToIndex(i-1, j, i ,j);
                m2=nodeToIndex(i, j, i,j+1);
                m3=nodeToIndex(i+1, j, i,j);
                m4=nodeToIndex(i, j-1, i ,j);
                
                if (m4 != 0){
                    A[index][m4-1] = 1.0;
                }if (m1 != 0){
                    A[index][m1-1] = 1.0;
                }if (m2 != 0){
                    A[index][m2-1] = -1.0;    
                }if (m3 != 0){
                    A[index][m3-1] = -1.0;
                }
                
                if ((i == iin) && (j == jin)){ //inicio
                    std::cout << "inicio: ";
                    b[index] = 1.0;
                    std::cout << " /inicio" <<std::endl;
                    ++index;
                }else if ((i == ifin) && (j == jfin)){ // final
                    std::cout << "final: ";

                    b[index] = -1.0;
                    std::cout << " /final" <<std::endl;
                    ++index;
                    
                }else if (!eliminado && isCorner(i, j, m, n)){
                    if (m4 != 0){
                    A[index][m4-1] = 0;
                    }if (m1 != 0){
                        A[index][m1-1] = 0;
                    }if (m2 != 0){
                        A[index][m2-1] = 0;    
                    }if (m3 != 0){
                        A[index][m3-1] = 0;
                    }
                    eliminado = true;
                }else{
                    ++index;
                }
            } // end for j
        }//end for i
        std::cout << "end nodos" << std::endl;
    }

    //ecuacioens de mallas 
    template <typename T>
    void mallas(int m, int n, Matrix<T> &A, const std::vector<int> &resistors){
        int index = m*n -1; //initial index to place the equations in A
        int m1, m2, m3, m4;
        for (int i = 0; i < (m-1); ++i){
            for (int j = 0; j< (n-1); ++j){
                //obtener posiciones para las resitencias
                m1 = nodeToIndex(i, j, i, j+1); // derecha adyacente al nodo (pos]
                m2 = nodeToIndex(i, j+1, i+1, j+1); //derecha para abajo (pos]
                m3 = nodeToIndex(i+1, j, i+1, j+1); // abajo acostada (neg]
                m4 = nodeToIndex(i , j, i+1, j); //abajo adyancente al nodo (neg]

                A[index][m1-1] = (T) resistors[m1-1];
                A[index][m2-1] = (T) resistors[m2-1];
                A[index][m3-1] = (T) -resistors[m3-1];
                A[index][m4-1] = (T) -resistors[m4-1];    
                index++;       
            }
        }
    }

    /**
     * Finds the trajectory of the biggest current 
     * rowI: initial row
     * colI: initial col 
     * */
    template <typename T>
    void trajectoryFinder(size_t rowI, size_t colI,size_t rowF, size_t colF, std::vector<size_t> x_){
        std::vector<size_t> trajectoryX;
        std::vector<size_t> trajectoryY;

        size_t n = rawMap_.rows();
        size_t m = rawMap_.cols();

        size_t colN = colI + 1; //next column
        size_t colP = colI - 1; //previous column
        size_t rowN = rowI + 1; //next row
        size_t rowP = rowI - 1; //previous row

        size_t checkedRow = rowI; //node already checked
        size_t checkedCol = colI;

        trajectoryX.push_back(rowI);
        trajectoryY.push_back(colI);

        size_t index1 = 0;
        size_t index2 = 0;
        size_t index3 = 0;

        size_t higher;

        indexPair pair;


        while(colI != colF && rowI != rowF){
            if(rowI < 0 || colI < 0){
                throw anpi::Exception("Punto inicial menor a cero");
            }
            else if((checkedRow == rowI && checkedCol == colI) && (colN != m && rowN != n) && colP > 0){
                index1 = nodeToIndex(rowI,colI,rowI,colN);
                index2 = nodeToIndex(rowI,colI,rowN,colI);
                index3 = nodeToIndex(rowI,colP,rowI,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,x_[index3],index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(checkedRow);
                trajectoryY.push_back(checkedCol);
            }
            else if((checkedRow == rowI && checkedCol == colI) && (colN != m && rowN != n) && rowP > 0){
                index1 = nodeToIndex(rowI,colI,rowI,colN);
                index2 = nodeToIndex(rowI,colI,rowN,colI);
                index3 = nodeToIndex(rowP,colI,rowI,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,x_[index3],index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(rowI);
                trajectoryY.push_back(colI);
            }
            else if((checkedRow == rowI && checkedCol == colI) && colN == m){
                index1 = nodeToIndex(rowI,colI,rowN,colI);
                index2 = nodeToIndex(rowI,colP,rowI,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(rowI);
                trajectoryY.push_back(colI);
            }
            else if((checkedRow == rowI && checkedCol == colI) && rowN == m){
                index1 = nodeToIndex(rowI,colI,rowI,colN);
                index2 = nodeToIndex(rowI,colP,rowI,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(rowI);
                trajectoryY.push_back(colI);
            }
            else if((checkedRow == rowI && checkedCol == colI) && colP < 0){
                index1 = nodeToIndex(rowI,colI,rowI,colN);
                index2 = nodeToIndex(rowI,colI,rowN,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(rowI);
                trajectoryY.push_back(colI);
            }
            else if((checkedRow == rowI && checkedCol == colI) && rowP < 0){
                index1 = nodeToIndex(rowI,colI,rowI,colN);
                index2 = nodeToIndex(rowI,colI,rowN,colI);
                higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                pair = indexToNodes(higher);
                rowI = pair.row1;
                colI = pair.col1;
                trajectoryX.push_back(rowI);
                trajectoryY.push_back(colI);
            }
            else{
                if((checkedRow == rowP && checkedCol == colP) && (colN != m && rowN != n) && colP > 0){
                    index1 = nodeToIndex(rowI,colI,rowI,colN);
                    index2 = nodeToIndex(rowI,colI,rowN,colI);
                    index3 = nodeToIndex(rowI,colP,rowI,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,x_[index3],index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
                else if((checkedRow == rowP && checkedCol == colP) && (colN != m && rowN != n) && rowP > 0){
                    index1 = nodeToIndex(rowI,colI,rowI,colN);
                    index2 = nodeToIndex(rowI,colI,rowN,colI);
                    index3 = nodeToIndex(rowP,colI,rowI,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,x_[index3],index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
                else if((checkedRow == rowP && checkedCol == colI) && colN == m){
                    index1 = nodeToIndex(rowI,colI,rowN,colI);
                    index2 = nodeToIndex(rowI,colP,rowI,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
                else if((checkedRow == rowP && checkedCol == colI) && rowN == m){
                    index1 = nodeToIndex(rowI,colI,rowI,colN);
                    index2 = nodeToIndex(rowI,colP,rowI,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
                else if((checkedRow == rowP && checkedCol == colI) && colP < 0){
                    index1 = nodeToIndex(rowI,colI,rowI,colN);
                    index2 = nodeToIndex(rowI,colI,rowN,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
                else if((checkedRow == rowP && checkedCol == colI) && rowP < 0){
                    index1 = nodeToIndex(rowI,colI,rowI,colN);
                    index2 = nodeToIndex(rowI,colI,rowN,colI);
                    higher = getHigher(x_[index1],index1,x_[index2],index2,0,index3);
                    pair = indexToNodes(higher);
                    checkedRow = rowI;
                    checkedCol = colI;
                    rowI = pair.row1;
                    colI = pair.col1;
                    trajectoryX.push_back(rowI);
                    trajectoryY.push_back(colI);
                }
            }

        }
    }

    /**
     * Give the higher number
     * */
    size_t getHigher(size_t a, size_t ia, size_t b, size_t ib, size_t c,size_t ic){
        size_t higher;
        if(a >= b){
            if(a >=c){
                higher = ia;
            }
            else{
                higher = ic;
            }
        }
        if(b >= a){
            if(b >=c){
                higher = ib;
            }
            else{
                higher = ic;
            }
        }
        if(c >= a){
            if(c >= b){
                higher = ic;
            }
            else{
                higher = ib;
            }
        }
        return higher;
    }


    /**
     * Construct the grid from the given file
     * @retturn true if successful or false otherwise 
     * */
    bool build(const std::string filename);

    /**
     * Compute the internal data to navigate between the given nodes
     * */
    bool navigate(const indexPair& nodes);
        
    };

} //anpi

#endif
