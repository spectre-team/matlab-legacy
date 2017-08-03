/*
* Dataset_opencvTest.cpp
* Tests Dataset_opencv
*
Copyright 2017 Wojciech Wilgierz

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


#include <gtest/gtest.h>
#include <span.h>
#include "Spectre.libClassifier/Dataset_opencv.h"
#include "Spectre.libException/InconsistentArgumentSizesException.h"
#include "Spectre.libException/OutOfRangeException.h"

namespace
{
using namespace Spectre::libClassifier;
using namespace Spectre::libException;

class Dataset_opencvInitializationTest : public ::testing::Test
{
public:
protected:
	const std::vector<DataType> data_long{ 0.5f, 0.4f, 0.6f, 1.1f, 1.6f, 0.7f, 2.1f, 1.0f, 0.6f };
	const std::vector<DataType> data_short{ 0.5f, 0.4f, 0.6f };
	const std::vector<Label> labels{ 3, 7, 14 };
	const std::vector<Label> labels_too_long{ 3, 7, 14, 5 };
};

TEST_F(Dataset_opencvInitializationTest, correct_dataset_opencv_initialization_col_size_one)
{
	EXPECT_NO_THROW(Dataset_opencv(data_short, labels));
}

TEST_F(Dataset_opencvInitializationTest, correct_dataset_opencv_initialization)
{
	EXPECT_NO_THROW(Dataset_opencv(data_long, labels));
}

TEST_F(Dataset_opencvInitializationTest, throws_for_inconsistent_size)
{
	EXPECT_THROW(Dataset_opencv(data_short, labels_too_long), InconsistentArgumentSizesException);
}

TEST_F(Dataset_opencvInitializationTest, correct_dataset_opencv_initialization_from_mat_col_size_one)
{
	std::vector<DataType> tmp_data(data_short);
	std::vector<Label> tmp_labels(labels);
	cv::Mat mat_data(3, 1, CV_TYPE, tmp_data.data());
	cv::Mat mat_labels(3, 1, CV_LABEL_TYPE, tmp_labels.data());
	EXPECT_NO_THROW(Dataset_opencv(mat_data, mat_labels));
}

TEST_F(Dataset_opencvInitializationTest, correct_dataset_opencv_initialization_from_mat)
{
	std::vector<DataType> tmp_data(data_long);
	std::vector<Label> tmp_labels(labels);
	cv::Mat mat_data(3, 3, CV_TYPE, tmp_data.data());
	cv::Mat mat_labels(3, 1, CV_LABEL_TYPE, tmp_labels.data());
	EXPECT_NO_THROW(Dataset_opencv(mat_data, mat_labels));
}

class Dataset_opencvTest : public ::testing::Test
{
public:
	Dataset_opencvTest()
	{
		
	}
protected:
	const std::vector<DataType> data{ 0.5f, 0.4f, 0.6f, 1.1f, 1.6f, 0.7f, 2.1f, 1.0f, 0.6f };
	const std::vector<Label> labels{ 3, 7, 14 };
	std::unique_ptr<Dataset_opencv> dataset;

	void SetUp() override
	{
		dataset = std::make_unique<Dataset_opencv>(data, labels);
	}
};

TEST_F(Dataset_opencvTest, get_out_of_range_data)
{
	EXPECT_THROW((*dataset)[4], OutOfRangeException);
}

TEST_F(Dataset_opencvTest, get_in_range_data)
{
	auto test = (*dataset)[1];
	std::vector<DataType> check({ 1.1f, 1.6f, 0.7f });
	EXPECT_EQ(test, check);
}

TEST_F(Dataset_opencvTest, get_out_of_range_label)
{
	EXPECT_THROW(dataset->GetSampleMetadata(4), OutOfRangeException);
}

TEST_F(Dataset_opencvTest, get_in_range_label)
{
	auto test = dataset->GetSampleMetadata(1);
	Label check = 7;
	EXPECT_EQ(test, check);
}

TEST_F(Dataset_opencvTest, get_dataset_metadata)
{
	Empty check = Empty::instance();
	EXPECT_EQ((*dataset).GetDatasetMetadata(), check);
}

TEST_F(Dataset_opencvTest, get_data)
{
	//nie potrafi porownac labels, wiec trzeba forem to robic
	EXPECT_EQ(dataset->GetData(), data);
}

TEST_F(Dataset_opencvTest, get_labels)
{
	//nie potrafi porownac labels, wiec trzeba forem to robic
	EXPECT_EQ(dataset->GetSampleMetadata(), labels);
}

TEST_F(Dataset_opencvTest, get_size)
{
	EXPECT_EQ(dataset->size(), labels.size());
}

TEST_F(Dataset_opencvTest, get_data_mat)
{
	std::vector<DataType> mat_data{ 0.5f, 0.4f, 0.6f, 1.1f, 1.6f, 0.7f, 2.1f, 1.0f, 0.6f };
	cv::Mat check(3, 3, CV_TYPE, mat_data.data());
	cv::Mat result = dataset->getMatData();
	EXPECT_EQ(check, result);
}

TEST_F(Dataset_opencvTest, get_labels_mat)
{
	std::vector<Label> mat_labels{ 3, 7, 14 };
	cv::Mat check(3, 1, CV_LABEL_TYPE, mat_labels.data());
	cv::Mat result = dataset->getMatLabels();
	EXPECT_EQ(check, result);
}

}
