// text i/o
#include <iostream>
#include <fstream>
// audio i/o
#include "AudioFile.h"
// image i/o
#include "EasyBMP.h"

// math
#include <cmath>
#include <complex>

// FFT
const double PI = acos(-1);

using namespace std;

template<typename T>
void fft(vector<complex<T>> & a, bool invert);

template<typename T>
void grab_audio_data(vector< vector<T> >& samples, int& num_channels, int& num_samples, int& sample_rate, const string filename);

template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, int frame);

void test_image();

void test_fft();

int main()
{
    vector< vector<double> > samples; int num_channels, num_samples, sample_rate;

    grab_audio_data(samples, num_channels, num_samples, sample_rate, "../a_lady.wav");
    double length_in_seconds = double(num_samples)/sample_rate;

    int frames_per_second = 12;
    double block_length_in_seconds = 1.0/frames_per_second;
    int block_length = block_length_in_seconds * sample_rate;
    cout << "Target FPS: " << frames_per_second << endl;
    cout << "Target block_length_in_seconds: " << block_length_in_seconds << endl;
    cout << "Target block_length: " << block_length << endl;
    block_length = pow(2, int(log2(block_length))); //rescale to an efficient power of 2
    block_length_in_seconds = block_length / double(sample_rate);
    frames_per_second = 1.0/block_length_in_seconds;
    cout << "Actual FPS: " << frames_per_second << endl;
    cout << "Actual block_length_in_seconds: " << block_length_in_seconds << endl;
    cout << "Actual block_length: " << block_length << endl;

    int num_blocks = int(num_samples / block_length); // How many images to generate
    cout << "Frame count: " << num_blocks << endl;

    vector<complex<double>> block_audio (block_length, 0);
    for(int block=0; block<num_blocks; block++)
    {
        for(int i=0; i<block_length; i++)
            block_audio[i] = 1.0*samples[0][block*block_length+i] + 1.0i*samples[1][block*block_length+i]; //left channel real, right channel imaginary....

        fft(block_audio, false);

        output_graph(block_audio, 1800, 900, block);

        cout << "frame " << block << " done!" << endl;
    }

    // Tests

    int current_sample;
    double stt_sec = 1.0;
    double end_sec = 1.01;
    for(int sample=int(stt_sec*sample_rate); sample<int(end_sec*sample_rate); sample+=5) // stt_sec - end_sec seconds of audio
    {
        cout << "|";
        for(int channel=0; channel<num_channels; channel++) // all channels
        {
            current_sample = 20*samples[channel][sample];
            for(int i=0; i<10+current_sample; i++)
                cout << "_";
            cout << "0";
            for(int i=0; i<10-current_sample; i++)
                cout << "_";
            cout << "|";
        }
        cout << double(sample)/sample_rate << endl;
    }

    test_image();

    test_fft();

    return 0;
}

template<typename T>
void grab_audio_data(vector< vector<T> >& samples, int& num_channels, int& num_samples, int& sample_rate, const string filename)
{
    // Create an AudioFile object
    AudioFile<T> audioFile;

    // Load an audio file
    bool loadedOK = audioFile.load (filename);
    assert (loadedOK);

    samples = audioFile.samples;
    //audioFile.samples is vector<vector<type>> of shape (num_channels, num_samples)
    //audioFile.samples[channel][sampleIndex] is <type>

    // Get some information about the loaded audio
    num_channels = audioFile.getNumChannels(); // number of channels
    num_samples = audioFile.getNumSamplesPerChannel(); // total number of samples
    sample_rate = audioFile.getSampleRate(); // samples per second
    cout << "audioFile.getLengthInSeconds() = " << audioFile.getLengthInSeconds() << " =?= " << double(num_samples)/sample_rate << " = num_samples/sample_rate" << endl;
    cout << "audioFile.getBitDepth() = " << audioFile.getBitDepth() << " <?= " << 32 << " = float bit depth" << endl;

    // or, just use this quick shortcut to print a summary to the console
    audioFile.printSummary();

    return;
}

template<typename T>
void fft(vector<complex<T>> & a, bool invert) {
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        complex<T> wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            complex<T> w(1);
            for (int j = 0; j < len / 2; j++) {
                complex<T> u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (complex<T> & x : a)
            x /= n;
    }
}

string fixed_len(int x, int len)
{
    const auto digits = int(log10(abs(x)));
    string str = "";
    for(int i=0; i<len-digits; i++)
        str +="0";
    str += to_string(x);

    //cout << x << " has " << digits << " digits so " << str << " has " << len << " digits." << endl;
    return str;
}
template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, int frame)
{
    // Parameters
    int bit_depth = 32;
    // Declare a new bitmap object
    BMP image;
    // Set size to 640 × 480
    image.SetSize(wdth, hght);
    // Set its color depth to 8-bits
    image.SetBitDepth(bit_depth);

    // Change colors
    for(int x=0; x<wdth; x++)
        for(int y=0; y<hght; y++)
            image.SetPixel(x, y, {0, 0, 0, 255});

    int N = values.size();

    double max_val = 0;
    for(int i=0; i<N; i+=1)
        if(abs(values[i])>max_val)
            max_val = abs(values[i]);

    int i, yr, yi;
    for(int x=0; x<wdth; x++)
    {
        i = double(N)*x/wdth;
        yr = int(0.5*hght + 0.5*hght*values[i].real()/max_val);
        yi = int(0.5*hght + 0.5*hght*values[i].imag()/max_val);

        image.SetPixel(x, yr, {0, 0, 255, 255}); //red for real
        image.SetPixel(x, yi, {255, 0, 0, 255}); //blue for imaginary
        //cout << "x=" << x << " yr=" << yr << " yi=" << yi << endl;
    }
    // Save to file
    string frame_str = "../frames/"+fixed_len(frame, 5)+"graph.bmp";
    //cout << "wrote " << frame_str << endl;
    image.WriteToFile(frame_str.c_str());
}

void test_image()
{
    // Parameters
    int wdth = 1000;
    int hght = 1000;
    int bit_depth = 32;
    // Declare a new bitmap object
    BMP image;
    // Set size to 640 × 480
    image.SetSize(wdth, hght);
    // Set its color depth to 8-bits
    image.SetBitDepth(bit_depth);
    // Change colors
    int r, g, b, a=255; // test with 8 bit depth
    for(int x=0; x<wdth; x++)
    {
        for(int y=0; y<hght; y++)
        {
            r = (255.0*x)/wdth; g = 0; b = (255.0*y)/hght;
            //cout << r << " ";
            image.SetPixel(x, y, {b, g, r, a});
        } //cout << endl;
    }
    // Save to file
    image.WriteToFile("../test.bmp");

    cout << "Image test done." << endl;

    return;
}

void test_fft()
{
    int N = pow(2, 10);

    vector<complex<double>> input (N, 0);
    vector<complex<double>> output (N, 0);

    double y;
    for(int i=0; i<N; i++)
    {
        y = 5*sin(50*2*PI*i/N) + 4*sin(40*2*PI*i/N) + 3*sin(30*2*PI*i/N) + 2*sin(20*2*PI*i/N) + 1*sin(10*2*PI*i/N);
        input[i] = y + 0.0i;

        output[i] = y + 0.0i;
    }

    fft(output, false);
    //fft(output, true); //ondoes the FFT

    // Plot results to check
    double in, input_max = 0;
    double out, output_max = 0;
    for(int i=0; i<N; i+=1)
    {
        in = abs(input[i]); //abs value
        if(in>input_max) input_max = in;
        out = abs(output[i]); //abs value
        if(out>output_max) output_max = out;
    }

    int current_sample;
    for(int i=0; i<N; i+=1)
    {
        cout << "|";

        current_sample = 10.*input[i].real()/input_max;
        for(int i=0; i<10+current_sample; i++)
            cout << "_";
        cout << "0";
        for(int i=0; i<10-current_sample; i++)
            cout << "_";

        cout << "|";

        current_sample = 10.*input[i].imag()/input_max;
        for(int i=0; i<10+current_sample; i++)
            cout << "_";
        cout << "0";
        for(int i=0; i<10-current_sample; i++)
            cout << "_";

        cout << "|";

        current_sample = 10.*output[i].real()/output_max;
        for(int i=0; i<10+current_sample; i++)
            cout << "_";
        cout << "0";
        for(int i=0; i<10-current_sample; i++)
            cout << "_";

        cout << "|";

        current_sample = 10.*output[i].imag()/output_max;
        for(int i=0; i<10+current_sample; i++)
            cout << "_";
        cout << "0";
        for(int i=0; i<10-current_sample; i++)
            cout << "_";

        cout << "| " << i << endl;
    }

    return;
}
