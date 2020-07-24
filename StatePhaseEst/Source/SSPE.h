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
            //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
            //transform = eigenSolver.eigenvectors() * (eigenSolver.eigenvalues().cwiseSqrt().asDiagonal());
        }

        Eigen::VectorXd mean;
        Eigen::MatrixXd transform;

        Eigen::VectorXd operator()() const
        {
            static std::mt19937 gen{ std::random_device{}() };
            static std::normal_distribution<> dist;

            //return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
            return Eigen::VectorXd();
        }
    };

    class SSPE {
    public:
        SSPE() : 
            fs(500)
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

        bool setFreqs(Array<float> newFreqs)
        {
            freqs = vector(newFreqs);
            return true;
        }

        const int getFs()
        {
            return fs;
        }

        bool setDesFreqIndex(int foi_lowcut, int foi_highcut)
        {
            for (int i = 0; i < freqs.size(); i++)
            {
                if (freqs(i) > foi_lowcut && freqs(i) < foi_lowcut)
                {
                    foiIndex = i;
                    return freqsSet = true;
                }
            }
            return freqsSet = false;
        }

        Array<std::complex<double>> evalBuffer(const double* rpHistory, int histSize)
        { 
            //vector prevState;
            Array<vector> allX = Array<vector>();
            allX.add(prevState); // Append prevState

            //Matrix prevCov;
            Array<Matrix> allP = Array<Matrix>();
            allP.add(prevCov); // Append prevCov

            for (int i = 0; i < histSize; i++)
            {
                vector newState;
                Matrix newStateCov;
                oneStepKalman(newState, newStateCov, allX[i], rpHistory[i], allP[i], phi, Q, sigmaObs, M);
                allX.add(newState);
                allP.add(newStateCov);
            }

            allX.remove(0); // remove prevState
            allP.remove(0); // remove prevCov
            prevState = allX.getLast();
            prevCov = allP.getLast();

            Array<std::complex<double>> complexArray = Array<std::complex<double>>();
            for (int i = 0; i < histSize; i++)
            {
                complexArray.add(std::complex<double>(allX[i](foiIndex), allX[i](foiIndex+1)));
            }

            return complexArray;
            //return;
        }

        float prctile(float percent, Array<float> arr)
        {
            int nArr = arr.size();
            int pctIndex = std::ceil(nArr * percent / 100);
            return arr[pctIndex];
        }

        void genKalmanParams(Matrix& phi, Matrix& Q, Matrix& M, vector freqs, vector ampVec, vector sigmaFreqs)
        {
            
            float fs = 30000; // get from channel info
            int nFreqs = freqs.size();
            Array<Matrix> rotMat = Array<Matrix>(); // R(w)
            Array<Matrix> varMat = Array<Matrix>(); // Q

            for (int i = 0; i < nFreqs; i++)
            {
                // Add rotation matrix
                double freq = freqs(i);
                Matrix mat = Matrix(2, 2);
                mat(0, 0) = cos(2 * M_PI * freq / fs);
                mat(0, 1) = -sin(2 * M_PI * freq / fs);
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

            phi = Matrix(nFreqs * 2, nFreqs * 2);
            Q = Matrix(nFreqs * 2, nFreqs * 2);

            // Add rotMat and varMat to phi and Q along the diagonal to improve speed of processing
            for (int i = 0; i < nFreqs; i++)
            {
                // Zeroes before mat
                for (int j = 0; j < i * 2; j++)
                {
                    phi(i, j) = 0;
                    phi(i + 1, j) = 0;
                    Q(i, j) = 0;
                    Q(i + 1, j) = 0;
                }

                // Put mat along diagonal
                phi(i * 2, i * 2) = rotMat[i](0, 0);
                phi(i * 2, i * 2 + 1) = rotMat[i](1, 0);
                phi(i * 2 + 1, i * 2) = rotMat[i](0, 1);
                phi(i * 2 + 1, i * 2 + 1) = rotMat[i](1, 1);

                Q(i * 2, i * 2) = rotMat[i](0, 0);
                Q(i * 2, i * 2 + 1) = rotMat[i](1, 0);
                Q(i * 2 + 1, i * 2) = rotMat[i](0, 1);
                Q(i * 2 + 1, i * 2 + 1) = rotMat[i](1, 1);

                // Pad zeroes after mat
                for (int j = i * 2 + 2; j < nFreqs * 2; j++)
                {
                    phi(i * 2, j) = 0;
                    phi(i * 2 + 1, j) = 0;
                    Q(i * 2, j) = 0;
                    Q(i * 2 + 1, j) = 0;
                }
            }
            
        }

        void fitModel(vector data, float fs)
        {
            if (~freqsSet)
            {
                // Freqs have not been set yet. need to set before moving on.
                jassertfalse;
                return;
            }

            vector freqEst = freqs / fs;
            int nSamples = data.size();
            int nFreqs = freqs.size();

            vector xstart = vector::Zero(2 * nFreqs);
            Matrix Pi;
            Pi.setIdentity(2 * freqs.size(), 2 * freqs.size());
            Matrix P = 0.001 * Pi;

            Matrix tempPhi;
            Matrix tempQ;
            Eigen::MatrixXf tempM;
            float tempSigmaObs;
            vector ampEst;
            vector allQ;

            // Create initial phi, Q and M based on params used last ( or estimates if first run)
            if (hasBeenFit())
            {
                Matrix tempPhi = phi;
                Matrix tempQ = Q;
                Eigen::MatrixXf tempM = M;
                float tempSigmaObs = sigmaObs;
                vector ampEst = ampVec;
                vector allQ = sigmaFreqs;
            }
            else
            {
                // Run AR model

                // Just hardcode for now?
            }
            
            genKalmanParams(tempPhi, tempQ, tempM, freqEst, ampEst, allQ); // gen initial kalman params

            int iter = 1;
            float errorVal = FLT_MAX;
            vector prevStateModel = prevState;
            Matrix prevCovModel = prevCov;

            while (iter < 400 && errorVal > 1e-3)
            {
                // run forward kalman filter
                Array<vector> allX = Array<vector>();
                allX.add(prevStateModel); // Append prevState

                Array<Matrix> allP = Array<Matrix>();
                allP.add(prevCovModel); // Append prevCov

                for (int i = 0; i < nSamples; i++)
                {
                    vector newState;
                    Matrix newStateCov;
                    oneStepKalman(newState, newStateCov, allX[i], data(i), allP[i], tempPhi, tempQ, tempSigmaObs, tempM);
                    allX.add(newState);
                    allP.add(newStateCov);
                }
                allX.remove(0); // remove prevState
                allP.remove(0); // remove prevCov

                // Run fixed interval smoother - Shumway, Stoffer, 82
                Array<vector> allXSmooth = Array<vector>();
                Array<Matrix> allPSmooth = Array<Matrix>();
                Array<Matrix> allJSmooth = Array<Matrix>();

                for (int i = 1; i < nSamples; i++)
                {
                    vector smoothState;
                    Matrix smoothCov;
                    Matrix smoothJ;
                    smoother(smoothState, smoothCov, smoothJ, allX[i], allX[i - 1], allP[i], allP[i - 1], tempPhi, tempQ);
                    allXSmooth.add(smoothState);
                    allPSmooth.add(smoothCov);
                    allJSmooth.add(smoothJ);
                }

                // need to estimate P_t_(t-1) for t = 1:N
                Matrix PTemp = tempPhi * allP.getLast() * tempPhi.transpose() + tempQ;
                Eigen::ArrayXXf PArrayForm = (tempM.transpose() * PTemp * tempM);
                Matrix PTransformedMatrix = 1 / (PArrayForm + tempSigmaObs);
                vector K = PTemp * tempM * PTransformedMatrix;
                Matrix Identity = Matrix().Identity(PTemp.rows(), PTemp.cols());
                Matrix P_N_N1 = (Identity - K * tempM.transpose()) * tempPhi * allP.getLast();


                Array<Matrix> allP_N_N1 = Array<Matrix>();
                allP_N_N1.add(P_N_N1);

                for (int i = nSamples - 1; i > 0; i--)
                {
                    Matrix P_N_N1 = allP[i] * allJSmooth[i - 1].transpose() + allJSmooth[i] * (allP_N_N1[i + 1] - tempPhi * allP[i]) * allJSmooth[i - 1];
                    allP_N_N1.insert(0, P_N_N1);
                }

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
                        B = B + allP_N_N1[i] + allXSmooth[i] * allXSmooth[i - 1];
                    }
                    C = C + allPSmooth[i] + allXSmooth[i] * allXSmooth[i];

                    // MS - Yikes... So you can only add scalar values to Arrays for some reason.
                    // MS - Also it doesn't let you treat the product of two Matrixes as a double, even if it will be. so temp is a 1x1 matrix and we just pull out the single number
                    //R = R+ M'* squeeze(newAllP(:,:,i)) * M + (y(i) - M'*newAllX(:,i)) *(y(i) - M'*newAllX(:,i))' ;
                    Eigen::ArrayXXf temp = tempM.transpose() * allPSmooth[i] * tempM;
                    Eigen::ArrayXXf temp1 = tempM.transpose() * allXSmooth[i].transpose();
                    Eigen::ArrayXXf temp2 = (tempM.transpose() * allXSmooth[i]).transpose();
                    tempSigmaObs = tempSigmaObs + temp(0) + ((data(i) - temp1) * (data(i) - temp2))(0);
                }

                tempSigmaObs = (1 / nSamples) * tempSigmaObs;

                vector oldFreq = freqEst * 1000 / (2 * M_PI);

                freqEst = vector::Zero(nFreqs);
                vector ampEst = vector::Zero(nFreqs);
                vector allQ = vector::Zero(nFreqs);

                for (int i = 0; i < nFreqs; i++)
                {

                    Matrix ATmp = extract2x2(B, i);
                    Matrix BTmp = extract2x2(B, i);
                    Matrix CTmp = extract2x2(B, i);

                    freqEst(i) = (atan((BTmp(2, 1) - BTmp(1, 2)) / BTmp.trace()));
                    ampEst(i) = sqrt(pow((BTmp(2, 1) - BTmp(1, 2)), 2) + pow(BTmp.trace(), 2)) / ATmp.trace();
                    allQ(i) = 1 / (2 * nSamples) * (CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace());

                    //// ALERT CHECK AMP EST POW FUNCTION. WE NEED TO DO .^2 vs other ones are just ^2
                    // Matlab functions below
                    //ampEst(numFreqs) = sqrt((B_tmp(2,1) - B_tmp(1,2))^2 + trace(B_tmp)^2)/trace(A_tmp);
                    //allQ(numFreqs) = 1/(2*length(y)) * (trace(C_tmp) - ampEst(numFreqs).^2 * trace(A_tmp));

                }

                genKalmanParams(tempPhi, tempQ, tempM, freqEst, ampEst, allQ);

                vector omega = freqEst * 1000 / (2 * M_PI);

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
            // }

            hasBeenUsed = true;
            
        }

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

           // return Matrix::Zero(2, 2);
        }


        void smoother(vector& smoothX, Matrix& smoothP, Matrix& smoothJ, vector x, vector x_prev, vector P, vector P_prev, Matrix phi, Matrix Q)
        {
            // Matlab code uses x and x_n, swapped to use x=x_prev and x_n=x. Ditto for P

            Matrix P_one = phi * P_prev * phi.transpose() + Q;

            smoothJ = P_prev * phi.transpose() * P_one.inverse();
            smoothX = x_prev + smoothJ * (x - phi * x_prev);
            smoothP = P_prev + smoothJ * (P - P_one) * smoothJ.transpose();
            
        }

        void oneStepKalman(vector& newState, Matrix& newStateCov, vector state, float measuredVal, Matrix stateCov, Matrix phi, Matrix Q, float sigmaObs, Matrix M)
        {
            vector stateEst = phi * state;
            Matrix phiEst = phi * stateCov * phi.transpose() + Q;

            Matrix temp = ((M.transpose() * phiEst * M).array() + sigmaObs).matrix();
            Matrix kalmanGain = (phiEst * M) / temp(0);


            temp = M.transpose() * stateEst;
            newState = stateEst + kalmanGain * (measuredVal - temp(0));
            newStateCov = phiEst - kalmanGain * M.transpose() * phiEst;
        }


    private:
        bool hasBeenUsed = false;
        bool freqsSet = false;

        const int fs;

        // SSPE decs
        Matrix phi;
        Matrix Q;
        Matrix M;
        float sigmaObs;
        vector ampVec;
        vector sigmaFreqs;

        vector freqs;

        vector prevState;
        Matrix prevCov;

        int foiIndex;


        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SSPE);
    };
}

#endif // SSPE_H_INCLUDED
