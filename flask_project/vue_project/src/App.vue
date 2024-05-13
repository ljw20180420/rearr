<script setup>
import { ref } from 'vue';
import uploadFile from './components/uploadFile.vue';
import runJob from './components/runJob.vue';
import setValues from './components/setValues.vue';
const pageStatus = ref({
  "target file": null,
  "pair file": null,
  "file without duplicates": null,
  "target spliter": null,
  "target spliter index": null,
  "pair spliter": null,
  "pair spliter index": null,
  "demultiplex file": null,
  "file of reads to align": null,
  "file of reference": null,
  "alignments": null,
  "minimal alignment score of target spliter": null,
  "minimal alignment score of pair spliter": null,
  "mismatching score": -6,
  "matching score for non-extension reference part": 4,
  "matching score for extension reference part": 2,
  "gap-extending penalty": -3,
  "gap-opening penalty": -9,
  "gap-extending penalty for unaligned reference end": 0,
  "gap-opening penalty for unaligned reference end": 0,
  "gap-extending penalty for unaligned query part": 0,
  "gap-opening penalty for unaligned query part": -5,
  "PAM1": "NGG",
  "PAM2": "NGG"
});
const display = ref(["target file", "pair file", "file without duplicates", "target spliter", "target spliter index", "pair spliter", "pair spliter index", "demultiplex file", "file of reads to align", "file of reference", "alignments"]);
</script>

<template>
  <div class="outer">
    <img src="./assets/images/projectLogic.png">
  </div>

  <div class="outer">
    <div v-for="(val, key) in pageStatus">
      <span v-if="display.includes(key)">
        {{ key }}: {{ val }}
      </span>
    </div>
  </div>

  <div class="outer">
    <uploadFile fileType="target file" v-model="pageStatus" />
    <uploadFile fileType="pair file" v-model="pageStatus" />
    <runJob url="removeDuplicates" fileType="file without duplicates" v-model="pageStatus" />
  </div>

  <div class="outer">
    <uploadFile fileType="target spliter" v-model="pageStatus" />
    <runJob url="buildSpliter/target spliter" fileType="target spliter index" v-model="pageStatus" />
    <uploadFile fileType="pair spliter" v-model="pageStatus" />
    <runJob url="buildSpliter/pair spliter" fileType="pair spliter index" v-model="pageStatus" />
  </div>

  <div class="outer">
    <uploadFile fileType="file without duplicates" v-model="pageStatus" />
    <setValues :params="['minimal alignment score of target spliter', 'minimal alignment score of pair spliter']" v-model="pageStatus" />
    <runJob url="demultiplex" fileType="demultiplex file" v-model="pageStatus" />
  </div>

  <div class="outer">
    <uploadFile fileType="file of reads to align" v-model="pageStatus" />
    <uploadFile fileType="file of reference" v-model="pageStatus" />
    <setValues :params="['mismatching score', 'matching score for non-extension reference part', 'matching score for extension reference part', 'gap-extending penalty', 'gap-opening penalty', 'gap-extending penalty for unaligned reference end', 'gap-opening penalty for unaligned reference end', 'gap-extending penalty for unaligned query part', 'gap-opening penalty for unaligned query part']" :selectParams="['PAM1', 'PAM2']" :selectValues="[['NGG', 'CCN'], ['NGG', 'CCN']]" v-model="pageStatus" />
    <runJob url="rearrange" fileType="alignments" v-model="pageStatus" />
  </div>

  <div class="outer">
    <h1>Monitor Progresses</h1>
    <a href="/flower" target="_blank">Celery Flower</a>
  </div>

  <div class="outer">
    <h1>Shiny Applications</h1>
    <a href="/shiny" target="_blank">Shiny Server</a>
  </div>
</template>

<style>
    .outer {
        outline: 0px solid #CCC;
        border: 10px solid #DDD;
        background-color: #999;
        overflow: auto;
    }
</style>
