/*

*/

#ifndef SSPE_H_INCLUDED
#define SSPE_H_INCLUDED

#define _USE_MATH_DEFINES

#include <BasicJuceHeader.h>
#include <Eigen\Dense>
#include <math.h>
#include <random>
#include <complex>
#include <chrono>

//using namespace StatePhaseEst;

namespace StatePhaseEst
{ 
    using vector = Eigen::VectorXf;
    using Matrix = Eigen::MatrixXf;

    enum STATE_PARAM
    {
        FREQS,
    };

    // generate variables to guess on quality of phase esimtation
    struct normal_random_variable
    {
        
        normal_random_variable(Eigen::MatrixXd const& covar)
            : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
        {}

        normal_random_variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
            : mean(mean)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
            transform = eigenSolver.eigenvectors() * (eigenSolver.eigenvalues().cwiseSqrt().asDiagonal());
        }

        Eigen::VectorXd mean;
        Eigen::MatrixXd transform;

        Eigen::VectorXd operator()() const
        {
            static std::mt19937 gen{ std::random_device{}() };
            static std::normal_distribution<> dist;

            return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
            //return Eigen::VectorXd();
        }
        
    };

    class SSPE {
    public:
        SSPE() 
        :   fs          (1000)
        ,   foiIndex    (0)
        ,   sigmaObs    (0)
        ,   stride      (40000/fs)
        {}

        ~SSPE() { }

        void reset()
        {
            hasBeenUsed = false;
            freqsSet = false;
        }

        bool hasBeenFit()
        {
            return hasBeenUsed;
        }

        // returns true if successful.
        bool setParams(int paramIndex, float value)
        {
            return true;
        }

        bool setHistSize(int histSize)
        {
            histSize = histSize;

            allPEB = Array<Matrix>();
            allXEB = Array<Matrix>();

            for (int i = 0; i <= histSize; i+=stride)
            {
                allXEB.add(Matrix());
                allPEB.add(Matrix());
            }

            return true;
        }

        bool setFreqs(Array<float> newFreqs)
        {
            freqs = vector(newFreqs.size());
            for (int i = 0; i < newFreqs.size(); i++)
            {
                freqs(i) = newFreqs[i];
            }
            //freqs = vector(newFreqs);
            return true;
        }

        const int getFs()
        {
            return fs;
        }

        bool setDesFreqIndex(float foi_lowcut, float foi_highcut)
        {
            int nFreqs = freqs.size();
         
            for (int i = 0; i < nFreqs; i++)
            {
                if (freqs(i) > foi_lowcut && freqs(i) < foi_highcut)
                {
                    foiIndex = i;

                    freqsSet = true;
                    return freqsSet;
                }
            }

            freqsSet = false;
            return freqsSet;
        }

        Array<std::complex<double>> evalBuffer(const double* rpHistory, int histSize)
        { 
            //std::cout << "histSize: " << histSize << std::endl;
            //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            Array<Matrix> allX = Array<Matrix>();
            allX.set(0,prevState); // Append prevState

            Array<Matrix> allP = Array<Matrix>();
            allP.set(0,prevCov); // Append prevCov

            for (int i = 1; i <= histSize/stride; i++)
            {
                Matrix newState;
                Matrix newStateCov;
                oneStepKalman(newState, newStateCov, allX[i-1], rpHistory[(i-1)*stride], allP[i-1], phi, Q, sigmaObs, M);
                allX.set(i,newState);
                allP.set(i,newStateCov);
            }
            //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

            prevState = allX.getLast();
            prevCov = allP.getLast();

            Array<std::complex<double>> complexArray = Array<std::complex<double>>();

            for (int i = 1; i <= histSize / stride; i++)
            {
                complexArray.add(std::complex<double>(allX[i](foiIndex, 0), allX[i](foiIndex+1,0)));
            }
        
            //Array<std::complex<double>> complexArray = Array<std::complex<double>>();
            //std::cout << "Time difference evalbuffer kalman = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
            //std::cin.get();
            return complexArray;      
        }

        float prctile(float percent, Array<float> arr)
        {
            
            int nArr = arr.size();
            int pctIndex = std::ceil(nArr * percent / 100);
            return arr[pctIndex];
            
            return 0;
        }

        void genKalmanParams(Matrix& phi, Matrix& Q, Matrix& M, vector freqs, vector ampVec, vector sigmaFreqs)
        {
            
            //float fs = 30000; // get from channel info
            int nFreqs = freqs.size();
            Array<Matrix> rotMat = Array<Matrix>(); // R(w)
            Array<Matrix> varMat = Array<Matrix>(); // Q

            for (int i = 0; i < nFreqs; i++)
            {
                // Add rotation matrix
                double freq = freqs(i);
                Matrix mat = Matrix(2, 2);
                mat(0, 0) = cos(2 * M_PI * freq / fs);
                mat(0, 1) = -1 * sin(2 * M_PI * freq / fs);
                mat(1, 0) = sin(2 * M_PI * freq / fs);
                mat(1, 1) = cos(2 * M_PI * freq / fs);
                Matrix ampMat = mat * ampVec(i);
                rotMat.add(ampMat);

                // Add variance matrix
                mat = Matrix(2, 2);
                mat(0, 0) = sigmaFreqs(i);
                mat(0, 1) = 0;
                mat(1, 0) = 0;
                mat(1, 1) = sigmaFreqs(i);
                varMat.add(mat);
            }

            phi = Matrix::Zero(nFreqs * 2, nFreqs * 2);
            Q = Matrix::Zero(nFreqs * 2, nFreqs * 2);

            // Add rotMat and varMat to phi and Q along the diagonal to improve speed of processing
            for (int i = 0; i < nFreqs; i++)
            {
                // Zeroes before mat
                for (int j = 0; j < i * 2; j++)
                {
                    //phi(i, j) = 0;
                    //phi(i + 1, j) = 0;
                    //Q(i, j) = 0;
                    //Q(i + 1, j) = 0;
                }

                // Put mat along diagonal
                phi(i * 2, i * 2) = rotMat[i](0, 0);
                phi(i * 2, i * 2 + 1) = rotMat[i](0, 1);
                phi(i * 2 + 1, i * 2) = rotMat[i](1, 0);
                phi(i * 2 + 1, i * 2 + 1) = rotMat[i](1, 1);

                Q(i * 2, i * 2) = varMat[i](0, 0);
                Q(i * 2, i * 2 + 1) = varMat[i](0, 1);
                Q(i * 2 + 1, i * 2) = varMat[i](1, 0);
                Q(i * 2 + 1, i * 2 + 1) = varMat[i](1, 1);

                // Pad zeroes after mat
                for (int j = i * 2 + 2; j < nFreqs * 2; j++)
                {
                    //phi(i * 2, j) = 0;
                    //phi(i * 2 + 1, j) = 0;
                    //Q(i * 2, j) = 0;
                    //Q(i * 2 + 1, j) = 0;
                }
            }
            
            //M = reshape([ones(length(freqs), 1), zeros(length(freqs), 1)]', length(freqs)*2,1); % for extracting real part of signal from analytic 'x'
            M = Matrix::Zero(nFreqs * 2,1);
            for (int i = 0; i < nFreqs * 2; i++)
            {
                if (i % 2 == 0)
                {
                    M(i,0) = 1;
                }
                //else
                //{
                    //M(i,0) = 0.0;
                //}
            }

        }

        void fitModel(Array<double> data_array)
        {
            std::cout << "Trying to fit" <<  std::endl;
            
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            int stride = 40000 / fs;
            const double* inputSeries = data_array.begin();
            vector data = vector(data_array.size() / stride);
            for (int i = 0; i < data_array.size()/stride; i++)
            {
                data(i) =  inputSeries[i*stride];
            }

            if (freqsSet == false)
            {
                // Freqs have not been set yet. need to set before moving on.
                jassertfalse;
                return;
            }
            
            int nSamples = data.size();
            int nFreqs = freqs.size();
            Matrix xstart = Matrix::Zero(2 * nFreqs, 1);
            Matrix Pi;
            Pi.setIdentity(2 * nFreqs, 2 * nFreqs);
            Matrix P = 0.001 * Pi;

            Matrix tempPhi;
            Matrix tempQ;
            Matrix tempM;
            float tempSigmaObs;
            vector ampEst;
            vector allQ;
            vector freqEst;
            vector omega;
            
            // Create initial phi, Q and M based on params used last ( or estimates if first run)
            if (hasBeenFit())
            {
                tempPhi = phi;
                tempQ = Q;
                tempM = M;
                tempSigmaObs = sigmaObs;
                ampEst = ampVec;
                allQ = sigmaFreqs;
                freqEst = freqs / fs;
            }
            else
            {
                // Run AR model

                // Just hardcode for now?
                //initParams(init_f, init_a, init_sigma, init_sigma_obs);
                /*
                ampEst = vector(1); // = // FROM MATLAB (0.9991)
                ampEst(0) = 0.999995434106301;
                allQ = vector(1); // = // from matlab (30.8456)
                allQ(0) = 20.8601143457368;
                tempSigmaObs = 0.0784202261595996; // = // from matlab
                freqEst = vector(1);
                freqEst(0) = 0.000628432315325676 / fs;
                freqs(0) = 0.000628432315325676;
                */

                ampEst = vector(3); // = // FROM MATLAB (0.9991)
                ampEst(0) = 0.99;
                ampEst(1) = 0.99;
                ampEst(2) = 0.99;
                allQ = vector(3); // = // from matlab (30.8456)
                allQ(0) = 50;
                allQ(1) = 50;
                allQ(2) = 50;
                tempSigmaObs = 1; // = // from matlab
                
                freqs = vector(3);
                freqs(0) = .3;
                freqs(1) = 10;
                freqs(2) = 60;
                
                freqEst = vector(3);
                freqEst = freqs / fs;
                nFreqs = freqs.size();
                Matrix Pi;
                Pi.setIdentity(2 * nFreqs, 2 * nFreqs);
                prevCov = 0.001 * Pi;
                P = 0.001 * Pi;
                //prevCov = Matrix::Zero(2*freqs.size(), 2*freqs.size()); // init prevCov and prevState to 0s
                prevState = Matrix::Zero(2 * nFreqs, 1);
                //tempQ = P;
                xstart = Matrix::Zero(2 * nFreqs, 1);
            }

            //vector freqEst = freqs / fs;
            
            genKalmanParams(tempPhi, tempQ, tempM, freqs, ampEst, allQ); // gen initial kalman params
            /*
            std::cout << "tempPhi:\n " << tempPhi << std::endl;
            std::cout << "tempQ:\n " << tempQ << std::endl;
            std::cout << "tempM:\n " << tempM << std::endl;
            */

            int iter = 1;
            float errorVal = FLT_MAX;
            Matrix prevStateModel = prevState;
            Matrix prevCovModel = prevCov;

            Array<Matrix> allP_N_N1 = Array<Matrix>();
            for (int i = 0; i < nSamples; i++)
            {
                allP_N_N1.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
            }

            Array<Matrix> allXSmooth = Array<Matrix>();
            Array<Matrix> allPSmooth = Array<Matrix>();
            Array<Matrix> allJSmooth = Array<Matrix>();
            for (int i = 0; i < nSamples; i++)
            {
                allXSmooth.add(Matrix::Zero(2 * nFreqs, 1));
                allPSmooth.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
                if (i != nSamples-1)
                {
                    allJSmooth.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
                }
            }

            while (iter < 400 && errorVal > 1e-3)
            {
                // run forward kalman filter
                Array<Matrix> allX = Array<Matrix>();
                allX.add(prevStateModel); // Append prevState

                Array<Matrix> allP = Array<Matrix>();
                allP.add(prevCovModel); // Append prevCov

                for (int i = 0; i < nSamples; i++)
                {
                    Matrix newState;
                    Matrix newStateCov;
                    oneStepKalman(newState, newStateCov, allX[i], data(i), allP[i], tempPhi, tempQ, tempSigmaObs, tempM);
                    allX.add(newState);
                    allP.add(newStateCov);

                }
                allX.remove(0); // remove prevState
                allP.remove(0); // remove prevCov

                // Run fixed interval smoother - Shumway, Stoffer, 82
                allXSmooth.set(nSamples-1, allX[nSamples-1]);
                allPSmooth.set(nSamples-1, allP[nSamples-1]);
               // allJSmooth = Array<Matrix>();
                //allXSmooth= allX;
                //allPSmooth = allP;
               // allJSmooth = allP;
                //allJSmooth.remove(nSamples - 1);

                for (int i = nSamples - 2; i >= 0; i--)
                {
                    Matrix smoothState;
                    Matrix smoothCov;
                    Matrix smoothJ;
                    int nxs = allX.size();
                    int nps = allP.size();
                    smoother(smoothState, smoothCov, smoothJ, allXSmooth[i + 1], allX[i], allPSmooth[i + 1], allP[i], tempPhi, tempQ);
                    //allXSmooth[i] = smoothState;
                    //allPSmooth[i] = smoothCov;
                    //allJSmooth[i] = smoothJ;
                    allXSmooth.set(i, smoothState);
                    allPSmooth.set(i, smoothCov);
                    allJSmooth.set(i, smoothJ);

                }

                /*
                std::cout << "allXSmooth: \n" << allXSmooth[0] << std::endl;
                std::cout << "allXSmooth end: \n" << allXSmooth[nSamples - 1] << std::endl;
                std::cout << "allPSmooth: \n" << allPSmooth[0] << std::endl;
                std::cout << "allPSmooth end: \n" << allPSmooth[nSamples-1] << std::endl;
                std::cout << "allJSmooth: \n" << allJSmooth[0] << std::endl;
                std::cout << "allJSmooth end: \n" << allJSmooth[nSamples - 2] << std::endl;
                std::cin.get();
                */
                // need to estimate P_t_(t-1) for t = 1:N
                Matrix PTemp = tempPhi * allP[nSamples-1] * tempPhi.transpose() + tempQ;
                double PArrayForm = (tempM.transpose() * PTemp * tempM)(0);
                double PTransformedMatrix = 1.0 / (PArrayForm + tempSigmaObs);
                Matrix K = PTemp * tempM * PTransformedMatrix;
                Matrix Identity = Matrix().Identity(PTemp.rows(), PTemp.cols());
                Matrix P_N_N1 = (Identity - K * tempM.transpose()) * tempPhi * allP[nSamples - 1];
                
                allP_N_N1.set(nSamples-1, P_N_N1);

                for (int i = nSamples - 2; i > 0; i--)
                {
                    Matrix P_N_N1 = allP[i] * allJSmooth[i - 1].transpose() + allJSmooth[i] * (allP_N_N1[i + 1] - tempPhi * allP[i]) * allJSmooth[i - 1].transpose();
                    allP_N_N1.set(i, P_N_N1);
                }

                //std::cout << "allpnn1: " << allP_N_N1[0] << std::endl;
                //std::cout << "allpnn1 end: " << allP_N_N1[nSamples-1] << std::endl;


                // EM Alg from Soulat, 2019
                Matrix A = P + xstart * xstart.transpose();
                Matrix B = Matrix::Zero(allPSmooth[0].rows(), allPSmooth[0].cols());
                Matrix C = Matrix::Zero(allPSmooth[0].rows(), allPSmooth[0].cols());
                double tempSigmaObs = 0;

                for (int i = 0; i < nSamples; i++)
                {
                    if (i > 1)
                    {
                        A = A + allPSmooth[i - 1] + allXSmooth[i - 1] * allXSmooth[i - 1].transpose();
                        //std::cout << "all_PNNN1 cols: " << allP_N_N1[i].cols() << "\nall_PNN1 rows: " << allP_N_N1[i].rows() << std::endl;
                        //std::cout << "allx cols: " << allXSmooth[i].cols() << "\nallx rows: " << allXSmooth[i].rows() << std::endl;
                        B = B + allP_N_N1[i] + allXSmooth[i] * allXSmooth[i - 1].transpose();
                    }
                    C = C + allPSmooth[i] + allXSmooth[i] * allXSmooth[i].transpose();

                    // MS - Yikes... So you can only add scalar values to Arrays for some reason.
                    // MS - Also it doesn't let you treat the product of two Matrixes as a double, even if it will be. so temp is a 1x1 matrix and we just pull out the single number
                    //R = R+ M'* squeeze(newAllP(:,:,i)) * M + (y(i) - M'*newAllX(:,i)) *(y(i) - M'*newAllX(:,i))' ;
                    double temp = (tempM.transpose() * allPSmooth[i] * tempM)(0);
                    double temp1 = (tempM.transpose() * allXSmooth[i])(0);
                    double temp2 = (data(i) - temp1);
                    //Eigen::ArrayXXf temp2 = (tempM.transpose() * allXSmooth[i]).transpose();
                    tempSigmaObs = tempSigmaObs + temp + (data(i) - temp1) * temp2;
                }

                tempSigmaObs = (1 / nSamples) * tempSigmaObs;

                vector oldFreq = freqEst * 1000 / (2 * M_PI);

                freqEst = vector::Zero(nFreqs);
                vector ampEst = vector::Zero(nFreqs);
                vector allQ = vector::Zero(nFreqs);

                for (int i = 0; i < nFreqs; i++)
                {

                    Matrix ATmp = extract2x2(A, i);
                    Matrix BTmp = extract2x2(B, i);
                    Matrix CTmp = extract2x2(C, i);
                    /*
                    std::cout << "after genK ATmp at end:\n " << ATmp << std::endl;
                    std::cout << "after genK BTmp at end:\n " << BTmp << std::endl;
                    std::cout << "(BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace()):\n " << (BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace() << std::endl;
                    std::cout << "after genK CTmp at end:\n " << CTmp << std::endl;
                    */

                    freqEst(i) = (atan((BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace()));
                    ampEst(i) = sqrt(pow((BTmp(1, 0) - BTmp(0, 1)), 2) + pow(BTmp.trace(), 2)) / ATmp.trace();
                    allQ(i) = (1.0 / float(2.0 * nSamples)) * (CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace());

                    //// ALERT CHECK AMP EST POW FUNCTION. WE NEED TO DO .^2 vs other ones are just ^2
                    // Matlab functions below
                    //ampEst(numFreqs) = sqrt((B_tmp(2,1) - B_tmp(1,2))^2 + trace(B_tmp)^2)/trace(A_tmp);
                    //allQ(numFreqs) = 1/(2*length(y)) * (trace(C_tmp) - ampEst(numFreqs).^2 * trace(A_tmp));

                }

                omega = freqEst * 1000 / (2 * M_PI);

                genKalmanParams(tempPhi, tempQ, tempM, omega, ampEst, allQ);
                /*
                std::cout << "after genK allQ at end:\n " << allQ << std::endl;
                std::cout << "after genK tempPhi at end:\n " << tempPhi << std::endl;
                std::cout << "after genK ampEst at end:\n " << ampEst << std::endl;
                std::cout << "after genK freqEst at end:\n " << freqEst << std::endl;
                std::cout << "after genK tempQ at end:\n " << tempQ << std::endl;
                std::cout << "after genK omega at end:\n " << omega << std::endl;
                std::cin.get();
                */
                // Reset state vector
                allX.clearQuick();
                allX.add(allXSmooth.getLast()); // Append prevState

                allP.clearQuick();
                allP.add(allPSmooth.getLast()); // Append prevCov

                iter += 1;
                //  MAKE SURE cwiseAbs does what we think it does
                errorVal = ((omega - oldFreq).cwiseAbs()).sum();
            }

            // Got good values send them to main thread. Maybe do this in run function

            // getLock()
            // {
            phi = tempPhi;
            Q = tempQ;
            M = tempM;
            sigmaObs = tempSigmaObs;
            ampVec = ampEst;
            sigmaFreqs = allQ;
            freqs = omega;
            hasBeenUsed = true;

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();


            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "M:\n" << M << std::endl;
            std::cout << "sigmaObs:\n" << sigmaObs << std::endl;
            std::cout << "ampVev:\n" << ampVec << std::endl;
            std::cout << "sigmaFreqs:\n" << sigmaFreqs << std::endl;
            std::cout << "freqs:\n" << freqs << std::endl;
            std::cout << "iter: " << iter << std::endl;
            
            // }

            
            std::cout << "Time difference fit model = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
            std::cin.get();
            
        }

        //Eigen::MatrixXf extract2x2(Matrix inputMatrix, int i)
        Eigen::MatrixXf extract2x2(Matrix inputMatrix, int i)
        {
            
            int rowcol = i * 2;
            Matrix x = Matrix::Zero(2, 2);
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    x(j, k) = inputMatrix(rowcol + j, rowcol + k);
                }
            }

            return x;
            
            //return Matrix::Zero(2, 2);
        }


        void smoother(Matrix& smoothX, Matrix& smoothP, Matrix& smoothJ, Matrix x, Matrix x_prev, Matrix P, Matrix P_prev, Matrix phi, Matrix Q)
        {
            // Matlab code uses x and x_n, swapped to use x=x_prev and x_n=x. Ditto for P
            Matrix P_one = phi * P_prev * phi.transpose() + Q;

            smoothJ = P_prev * phi.transpose() * P_one.inverse();
            smoothX = x_prev + smoothJ * (x - phi * x_prev);
            smoothP = P_prev + smoothJ * (P - P_one) * smoothJ.transpose();

            /*
            std::cout << "smoothX:\n" << smoothX << std::endl;
            std::cout << "smoothP:\n" << smoothP << std::endl;
            std::cout << "smoothJ:\n" << smoothJ << std::endl;
            std::cout << "x:\n" << x << std::endl;
            std::cout << "x_prev:\n" << x_prev << std::endl;
            std::cout << "P:\n" << P << std::endl;
            std::cout << "P_prev:\n" << P_prev << std::endl;
            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "P_one:\n" << P_one << std::endl;
            std::cin.get(); 
            */
        }

        void oneStepKalman(Matrix& newState, Matrix& newStateCov, Matrix state, float measuredVal, Matrix stateCov, Matrix phi, Matrix Q, float sigmaObs, Matrix M)
        {
            //float breakme = phi(13, 5364);
            Matrix stateEst = phi * state;
            Matrix pEst = phi * stateCov * phi.transpose() + Q;

            //int ncol = M.size();
            //int nlegn = phiEst.size();

            //float temp = (M.transpose() * pEst * M)(0) + sigmaObs;
            Matrix kalmanGain = (pEst * M) / ((M.transpose() * pEst * M)(0) + sigmaObs);

            Matrix temp = M.transpose() * stateEst;
            newState = stateEst + kalmanGain * (measuredVal - temp(0));
            newStateCov = pEst - kalmanGain * M.transpose() * pEst;

            /*
            std::cout << "newState:\n" << newState << std::endl;
            std::cout << "newStateCov:\n" << newStateCov << std::endl;
            std::cout << "state:\n" << state << std::endl;
            std::cout << "measuredVal:\n" << measuredVal << std::endl;
            std::cout << "stateCov:\n" << stateCov << std::endl;
            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "sigmaObs:\n" << sigmaObs << std::endl;
            std::cout << "M:\n" << M << std::endl;
            std::cout << "kalmanGain:\n" << kalmanGain << std::endl;
            std::cout << "stateEst:\n" << stateEst << std::endl;
            std::cout << "pEst:\n" << pEst << std::endl;
            std::cin.get();
            */
        }


    private:
        bool hasBeenUsed = false;
        bool freqsSet = false;

        const int fs;
        int histSize;
        int stride;

        // SSPE decs
        
        Matrix phi;
        Matrix Q;
        vector M;
       
        Array<Matrix> allXEB;

        Array<Matrix> allPEB;
        
        vector ampVec;
        vector sigmaFreqs;

        vector freqs;

        Matrix prevState;
        Matrix prevCov;
         
        int foiIndex;
        float sigmaObs;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SSPE);
    };
}

#endif 
