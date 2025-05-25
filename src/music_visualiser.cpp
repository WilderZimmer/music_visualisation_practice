// text i/o
#include <iostream>
// audio i/o
#include "AudioFile.h"
// image i/o
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "EasyBMP.h"
//command line tools
#include <cstdlib>

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

string fixed_len(int x, int len) // only positive x, lol !
{
    const auto digits = int(log10(abs(x)));
    string str = "";
    for(int i=0; i<len-1-digits; i++) str +="0";
    str += to_string(x);
    //cout << x << " has " << digits << " digits so " << str << " has " << len << " digits." << endl;
    return str;
}
template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, const string& frame_string);

void test_image();

void test_fft();

int main()
{
    test_image();

    string audio = "..\\wav_inputs\\a_lady.wav";
    string video = "..\\outputs\\a_lady.mp4";
    string video_with_sound = "..\\outputs\\a_lady_with_sound.mp4";

    vector< vector<double> > samples; 
    int num_channels, num_samples, sample_rate;
    grab_audio_data(samples, num_channels, num_samples, sample_rate, audio);
    double length_in_seconds = double(num_samples)/sample_rate;

    double frames_per_second = 12.0;
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

    //TODO buffer the calculation int one big vec<vec< cplx<dbl> >> with a ring buffer
    //TODO separate frame writing (slow) from calculation (fast)
    vector<complex<double>> block_audio (block_length, 0); 
    int str_len = log10(abs(num_blocks)) + 1;
    for(int block=0; block<num_blocks; block++)
    {
        for(int i=0; i<block_length; i++)
            block_audio[i] = 1.0*samples[0][block*block_length+i] + 1.0i*samples[1][block*block_length+i]; //left channel real, right channel imaginary....

        fft(block_audio, false);

        output_graph(block_audio, 1800, 900, "..\\output_frames\\"+fixed_len(block, str_len)+"graph.bmp");

        cout << "frame " << block << " done!" << endl;
    }

    cout << "Attempting ffmpeg..." << endl;

    string temp = "";
    int ffmpeg_fail_code;
    //Compile Frames with sequantial filenames <PRE><NUM><POST> (where NUM is a number of width N):
    temp = "ffmpeg -y -framerate 21.5 -i ..\\output_frames\\%0" + to_string(str_len) + "d" + "graph.bmp " + video + " 2>&1"; 
    ffmpeg_fail_code = system(temp.c_str()); 
    if (ffmpeg_fail_code) cerr << "ffmpeg video creation failed with code " << ffmpeg_fail_code << endl;
    //Add audio:
    temp = "ffmpeg -y -i " + video + " -i " + audio + " -map 0:v -map 1:a -c:v copy -shortest " + video_with_sound + " 2>&1";
    ffmpeg_fail_code = system(temp.c_str()); 
    if (ffmpeg_fail_code) cerr << "ffmpeg audio adding failed with code " << ffmpeg_fail_code << endl;

    cout << flush; cin.clear(); // reset in case of weird console state from system calls
    cout << "ffmpeg conversion done!" << endl;

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

template<typename T>
void output_graph(vector<complex<T>> & values, int wdth, int hght, const string& frame_str)
{
    // Parameters
    int bit_depth = 32;
    // Declare a new bitmap object
    BMP image;
    // Set size to 640 � 480
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
    //if(max_val == 0.0) max_val = 1.0;

    int i, yr, yi;
    for(int x=0; x<wdth; x++)
    {
        if (max_val > 0) // most cases
        {
            i = 0.5*double(N)*x/wdth; //transform from image scale to vector index (halved because symmetric)
            yr = int(0.5*hght + 0.5*hght*values[i].real()/max_val);
            yi = int(0.5*hght + 0.5*hght*values[i].imag()/max_val);
            if (abs(yr)>=hght) yr=int(((yr>0)-(yr<0)+1)*hght/2.0); // If I've fucked up and overflowed somewhere
            if (abs(yi)>=hght) yi=int(((yi>0)-(yi<0)+1)*hght/2.0); // If I've fucked up and overflowed somewhere
        }
        else // sometimes there's literally no sound
        {
            yr = 0.0;
            yi = 0.0;
        }

        image.SetPixel(x, yr, {0, 0, 255, 255}); //red for real
        image.SetPixel(x, yi, {255, 0, 0, 255}); //blue for imaginary

        for (int h=0; h<abs(yr-hght/2); h++)
            image.SetPixel(x, int(0.5*hght+((yr-hght/2>0)-(yr-hght/2<0))*h), {0, 0, 255, 255}); //draw red line
        for (int h=0; h<abs(yi-hght/2); h++) 
            image.SetPixel(x, int(0.5*hght+((yi-hght/2>0)-(yi-hght/2<0))*h), {255, 0, 0, 255}); //draw blue line
    }

    // Save to file
    //cout << "wrote " << frame_str << endl;
    image.WriteToFile(frame_str.c_str());
}

void set_buffer_pixel(vector<uint8_t> & image_buffer, int wdth, int x, int y, uint8_t r, uint8_t g, uint8_t b)
{
    //Manual indexing :(
    int idx = (y*wdth + x) * 3;
    image_buffer[idx + 0] = r;
    image_buffer[idx + 1] = g;
    image_buffer[idx + 2] = b;
    return;
}
void test_image()
{
    cout << "generating test image" << endl;
    // Parameters
    int wdth = 1024;
    int hght = 1024;
    // Declare a new buffer for a full image
    vector<uint8_t> image_buffer(wdth * hght * 3); // x * y * rgb // alternately do wdth*hght*3 for rgba
    // Change colors
    uint8_t r, g, b;
    for(int x=0; x<wdth; x++)
    {
        for(int y=0; y<hght; y++)
        {
            r = (255.0*x)/wdth; g = 0; b = (255.0*y)/hght;
            //cout << r << " ";
            set_buffer_pixel(image_buffer, wdth, x, y, r, g, b);
        } //cout << endl;
    }
    // Save to file
    stbi_write_png("..\\outputs\\test.png", wdth, hght, 3, image_buffer.data(), wdth * 3);

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

    cout << "testing fft..." << endl;
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
    cout << "fft test done!";

    return;
}
