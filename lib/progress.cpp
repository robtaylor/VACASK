#include "progress.h"
#include "common.h"


namespace NAMESPACE {

ProgressReporter::ProgressReporter(double dt) 
    : enabled_(true), extent_(1), pos_(0), dt(dt), prefix_(""), postfix_(""), 
      valueFormat_(ValueFormat::Default), precision_(6) {
}

}
