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
#include <string>
#include <cstdlib>

#include <Matrix.hpp>
#include <Exception.hpp>



#include <AnpiConfig.hpp>
#include <opencv2/core.hpp>    // For cv::Mat
#include <opencv2/highgui.hpp> // For cv::imread/imshow
#include <opencv2/imgcodecs/imgcodecs_c.h>
#include <opencv2/highgui/highgui_c.h>



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
    
    
    //get rawMap_ for build() tests
    Matrix<float> getRawMap(){
        return rawMap_;

    }
    ///Construct
    ResistorGrid(Matrix<float> A, std::vector<float> b, Matrix<float> rawMap){
        A_ = A;
        b_ = b;
        rawMap_ = rawMap;
    } 

    /**
     * Map the index of the terminal nodes
     * */
    std::size_t nodesToIndex(const std::size_t row1,
                            const std::size_t col1,
                            const std::size_t row2,
                            const std::size_t col2)
    {
        std::size_t result=-1;//variable para contener la respuesta
        if (row1 == row2 && col1 == col2)
        {
            throw anpi::Exception("Indices iguales");
        }
        else if (!((((row1+1) == row2 && col1 == col2) || ((row1-1) == row2 && col1 == col2)) || 
                ((row1 == row2 && (col1+1) == col2) || (row1 == row2 && (col1-1) == col2))))
        {
            throw anpi::Exception("Indices no son adyacentes");
        }
        else{
            if ( ((row1-1) == row2 && col1 == col2) || (row1 == row2 && (col1-1) == col2) )
            {
                throw anpi::Exception("Indices ingresados en orden incorrecto. Primero el de más arriba o más a la izquierda.");
            }
            else 
            {
               
                std::size_t m,n;
                m=rawMap_.cols();
                n=rawMap_.rows();
                if(row1 == row2)//nodos en la misma fila 
                {
                    result=row1*(m-1)+col1;
                }
                else if (col1 == col2)//nodos en la misma columna
                {
                    result=(n*m-n)+col1*(n-1)+row1;
                }
            }

        }
        return result;
    }   

    /// Convert an inde x to the pair of node coordinates
    indexPair indexToNodes (const std::size_t idx){
        anpi::indexPair coordenadas;
        size_t i1, i2, j1, j2;
        size_t m=rawMap_.cols();
        size_t n=rawMap_.rows();
        size_t rf=n*m-n;//referencia de ultimo indice utilizado en el mapeo horizontal

        if (idx <= rf)//indice corrresponde a una fila de resistencias
        {
            if ((idx+1)%(m-1)!=0)// definicion de fila del primer indice
            {
                i1=(idx+1)/(m-1);
            }else
            {
                i1=(idx+1)/(m-1)-1;
            }
            //definicion de la columna del primer indice
            j1=idx-i1*(m-1);
            //definicion de fila y columna del segundo indice
            i2=i1;
            j2=j1+1;
        }
        else//indice correspondiante a una columna de resistencias
        {
            if((idx-rf)%(n-1)!=0)//definicion de columna del primer indice
            {
                j1=(idx-rf)/(n-1);
            }else 
            {
                j1=(idx+rf)/(n-1)-1;
            }
            //definicion de fila del primer indice
            i1=(idx-rf)-j1*(n-1);
            //definicion de fila y columna del segundo indice
            j2=j1;
            i2=i1+1;
        }
        coordenadas={i1,j1,i2,j2};

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
                
                m1=nodesToIndex(i-1, j, i ,j);
                m2=nodesToIndex(i, j, i,j+1);
                m3=nodesToIndex(i+1, j, i,j);
                m4=nodesToIndex(i, j-1, i ,j);
                
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
                m1 = nodesToIndex(i, j, i, j+1); // derecha adyacente al nodo (pos)
                m2 = nodesToIndex(i, j+1, i+1, j+1); //derecha para abajo (pos)
                m3 = nodesToIndex(i+1, j, i+1, j+1); // abajo acostada (neg)
                m4 = nodesToIndex(i , j, i+1, j); //abajo adyancente al nodo (neg)

                A[index][m1-1] = (T) resistors[m1-1];
                A[index][m2-1] = (T) resistors[m2-1];
                A[index][m3-1] = (T) -resistors[m3-1];
                A[index][m4-1] = (T) -resistors[m4-1];    
                index++;       
            }
        }
    }

    /**
     * Construct the grid from the given file
     * @retturn true if successful or false otherwise 
     * */
    bool build(const std::string filename)
    {
        try{
        std::string mapPath = std::string(ANPI_DATA_PATH) + filename;
        // Read the image using the OpenCV
        cv::Mat_<float> map;

        cv::imread(mapPath.c_str(),
                   CV_LOAD_IMAGE_GRAYSCALE)
            .convertTo(map, CV_32FC1);
        map /= 255.0f; // normalize image range to 0 .. 255

        // And create a window to show the image
        cv::namedWindow(mapPath, CV_WINDOW_NORMAL | CV_GUI_EXPANDED);
        cv::imshow(mapPath, map);

        // Convert the OpenCV matrix into an anpi matrix
        // We have to use the std::allocator to avoid an exact stride
        anpi::Matrix<float, std::allocator<float>> amapTmp(map.rows,
                                                           map.cols,
                                                           map.ptr<float>());

        // And transform it to a SIMD-enabled matrix
        anpi::Matrix<float> amap(amapTmp);
        rawMap_=amap;
        cv::waitKey();
        return true;
        }
        catch (cv::Exception e){
            std::cout<< "Error al cargar imagen"<<std::endl;
        }
        return false;
    }

    /**
     * Compute the internal data to navigate between the given nodes
     * */
    bool navigate(const indexPair& nodes);
    };

} //anpi

#endif
