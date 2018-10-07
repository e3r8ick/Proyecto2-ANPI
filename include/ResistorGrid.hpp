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
        ///Row mof the first node
        std::size_t row1;
        //Column of the first node 
        std::size_t col1;
        ///Row mof the second node
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
                ((row1 == row2 && (col1+1) == col2) || (row1 == row2 && (col1-1) == col2))))
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
