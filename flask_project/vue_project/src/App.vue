<script setup>
import { markRaw } from 'vue'
import { VueFlow } from '@vue-flow/core'
import { ref } from 'vue';
import dataTank from './components/dataTank.vue';
import runJob from './components/runJob.vue';

const nodeTypes = {
  dataTankNode: markRaw(dataTank),
  runJobNode: markRaw(runJob)
}

const nodes = ref([
  { id: 'target file', type: 'dataTankNode', position: { x: 950, y: 50 }, data: [ { type: 'file', taskId: null, name: 'target file', value: [null] } ] },
  { id: 'pair file', type: 'dataTankNode', position: { x: 950, y: 150 }, data: [ { type: 'file', taskId: null, name: 'pair file', value: [null] } ] },
  { id: 'removeDuplicates', type: 'runJobNode', position: { x: 1600, y: 100 }, data: {
    values: { 'target file': [null], 'pair file': [null] },
    taskIds: { 'target file': null, 'pair file': null },
  }},
  { id: 'file without duplicates', type: 'dataTankNode', position: { x: 1850, y: 100 }, data: [ { type: 'file', taskId: null, name: 'file without duplicates', value: [null] } ] },

  { id: 'target spliter', type: 'dataTankNode', position: { x: 950, y: 250 }, data: [ { type: 'file', taskId: null, name: 'target spliter', value: [null] } ] },
  { id: 'buildSpliter/target spliter', type: 'runJobNode', position: { x: 1600, y: 250 }, data: {
    values: { 'target spliter': [null] },
    taskIds: { 'target spliter': null },
  }},
  { id: 'target spliter index', type: 'dataTankNode', position: { x: 1850, y: 250 }, data: [
    { type: 'value', name: 'minimal alignment score of target spliter', value: null },
    { type: 'files', taskId: null, name: 'target spliter index', value: new Array(6).fill(null) },
  ] },

  { id: 'pair spliter', type: 'dataTankNode', position: { x: 950, y: 450 }, data: [ { type: 'file', taskId: null, name: 'pair spliter', value: [null] } ] },
  { id: 'buildSpliter/pair spliter', type: 'runJobNode', position: { x: 1600, y: 450 }, data: {
    values: { 'pair spliter': [null] },
    taskIds: { 'pair spliter': null }
  }},
  { id: 'pair spliter index', type: 'dataTankNode', position: { x: 1850, y: 450 }, data: [
    { type: 'value', name: 'minimal alignment score of pair spliter', value: null },
    { type: 'files', taskId: null, name: 'pair spliter index', value: new Array(6).fill(null) },
  ] },

  { id: 'demultiplex', type: 'runJobNode', position: { x: 2500, y: 350 }, data: {
    values: {
      'file without duplicates': [null],
      'target spliter index': new Array(6).fill(null),
      'pair spliter index': new Array(6).fill(null),
      'minimal alignment score of target spliter': null,
      'minimal alignment score of pair spliter': null,
    },
    taskIds: {
      'file without duplicates': null,
      'target spliter index': null,
      'pair spliter index': null,
    },
  } },
  { id: 'demultiplex file', type: 'dataTankNode', position: { x: 2750, y: 350 }, data: [{ type: 'file', taskId: null, name: 'demultiplex file', value: [null] }] },

  { id: "minimal base number after remove 5' spliter and 3' adapter from target", type: 'dataTankNode', position: { x:2750, y: 550 }, data: [{ type: 'value', name: "minimal base number after remove 5' spliter and 3' adapter from target", value: 30 }] },
  { id: 'sxPostProcess', type: 'runJobNode', position: { x: 3400, y: 450 }, data: {
    values: {
      'demultiplex file': [null],
      "minimal base number after remove 5' spliter and 3' adapter from target": 30,
    },
    taskIds: {
      'demultiplex file': null,
    },
  } },
  { id: 'file of reads to align', type: 'dataTankNode', position: { x: 4550, y: 650 }, data: [
    { type: 'value', name: 'gap-extending penalty for unaligned query part', value: 0 },
    { type: 'value', name: 'gap-opening penalty for unaligned query part', value: -5 },
    { type: 'file', taskId: null, name: 'file of reads to align', value: [null] },
  ] },

  { id: 'file of reference', type: 'dataTankNode', position: { x: 4550, y: 950 }, data: [
    { type: 'value', name: 'gap-extending penalty for unaligned reference end', value: 0 },
    { type: 'value', name: 'gap-opening penalty for unaligned reference end', value: 0 },
    { type: 'select', name: 'PAM1', value: 'NGG', options: ['NGG', 'CCN'] },
    { type: 'select', name: 'PAM2', value: 'NGG', options: ['NGG', 'CCN'] },
    { type: 'file', taskId: null, name: 'file of reference', value: [null] },
  ] },
  { id: 'align score settings', type: 'dataTankNode', position: {x: 4550, y: 350 }, data: [
    { type: 'value', name: 'mismatching score', value: -6 },
    { type: 'value', name: 'matching score for non-extension reference part', value: 4 },
    { type: 'value', name: 'matching score for extension reference part', value: 2 },
    { type: 'value', name: 'gap-extending penalty', value: -3 },
    { type: 'value', name: 'gap-opening penalty', value: -9 },
  ] },
  { id: 'rearrange', type: 'runJobNode', position: { x: 5200, y: 650 }, data: {
    values: {
      'file of reads to align': [null],
      'file of reference': [null],
      'gap-extending penalty for unaligned query part': null,
      'gap-opening penalty for unaligned query part': null,
      'gap-extending penalty for unaligned reference end': null,
      'gap-opening penalty for unaligned reference end': null,
      'PAM1': null,
      'PAM2': null,
      'mismatching score': null,
      'matching score for non-extension reference part': null,
      'matching score for extension reference part': null,
      'gap-extending penalty': null,
      'gap-opening penalty': null,
    },
    taskIds: {
      'file of reads to align': null,
      'file of reference': null,
    },
  } },
  { id: 'alignments', type: 'dataTankNode', position: { x: 5450, y: 650 }, data: [{ type: 'file', taskId: null, name: 'alignments', value: [null] }] },

  { id: 'genome', type: 'dataTankNode', position: { x: 2750, y: 1050 }, data: [{ type: 'file', taskId: null, name: 'genome', value: ["../genome/genome.fa"] }] },
  { id: 'indexGenome', type: 'runJobNode', position: { x: 3400, y: 1150 }, data: {
    values: { 'genome': [null] },
    taskIds: { 'genome': null },
  } },
  { id: 'genome index', type: 'dataTankNode', position: { x: 3650, y: 1150 }, data: [{ type: 'files', taskId: null, name: 'genome index', value: [
    '../genome/genome.1.bt2',
    '../genome/genome.2.bt2',
    '../genome/genome.3.bt2',
    '../genome/genome.4.bt2',
    '../genome/genome.rev.1.bt2',
    '../genome/genome.rev.2.bt2',
  ] }] },

  { id: 'csvfile', type: 'dataTankNode', position: { x: 50, y: 650 }, data: [{ type: 'file', taskId: null, name: 'csvfile', value: [null] }] },
  { id: 'extension length', type: 'dataTankNode', position: { x: 3650, y: 750 }, data: [
    { type: 'value', name: 'cleavage 1 extend upstream', value: 50 },
    { type: 'value', name: 'cleavage 1 extend downstream', value: 0 },
    { type: 'value', name: 'cleavage 2 extend upstream', value: 10 },
    { type: 'value', name: 'cleavage 2 extend downstream', value: 100 },
  ] },
  { id: 'getReference', type: 'runJobNode', position: { x: 4300, y: 1050 }, data: {
    values: {
      'csvfile': [null],
      'genome': [null],
      'genome index': [null],
      'cleavage 1 extend upstream': null,
      'cleavage 1 extend downstream': null,
      'cleavage 2 extend upstream': null,
      'cleavage 2 extend downstream': null,
    },
    taskIds: {
      'csvfile': null,
      'genome': null,
      'genome index': null,
    },
  } },

  { id: 'getSpliters', type: 'runJobNode', position: {x: 700, y: 350}, data: {
    values: { 'csvfile': [null] },
    taskIds: { 'csvfile': null },
  } }
])

const edges = ref([
  { id: 'e0', source: 'target file', target: 'removeDuplicates' },
  { id: 'e1', source: 'pair file', target: 'removeDuplicates' },
  { id: 'e2', source: 'removeDuplicates', target: 'file without duplicates' },

  { id: 'e3', source: 'target spliter', target: 'buildSpliter/target spliter' },
  { id: 'e4', source: 'buildSpliter/target spliter', target: 'target spliter index' },

  { id: 'e5', source: 'pair spliter', target: 'buildSpliter/pair spliter' },
  { id: 'e6', source: 'buildSpliter/pair spliter', target: 'pair spliter index' },

  { id: 'e7', source: 'file without duplicates', target: 'demultiplex' },
  { id: 'e8', source: 'target spliter index', target: 'demultiplex' },
  { id: 'e9', source: 'pair spliter index', target: 'demultiplex' },
  { id: 'e10', source: 'demultiplex', target: 'demultiplex file' },

  { id: 'e11', source: 'demultiplex file', target: 'sxPostProcess' },
  { id: 'e12', source: "minimal base number after remove 5' spliter and 3' adapter from target", target: 'sxPostProcess' },
  { id: 'e13', source: 'sxPostProcess', target: 'file of reads to align' },

  { id: 'e14', source: 'file of reads to align', target: 'rearrange' },
  { id: 'e15', source: 'file of reference', target: 'rearrange' },
  { id: 'e16', source: 'align score settings', target: 'rearrange' },
  { id: 'e17', source: 'rearrange', target: 'alignments' },

  { id: 'e18', source: 'genome', target: 'indexGenome' },
  { id: 'e19', source: 'indexGenome', target: 'genome index' },

  { id: 'e20', source: 'csvfile', target: 'getReference' },
  { id: 'e21', source: 'genome', target: 'getReference' },
  { id: 'e22', source: 'genome index', target: 'getReference' },
  { id: 'e23', source: 'extension length', target: 'getReference' },
  { id: 'e24', source: 'getReference', target: 'file of reference' },

  { id: 'e25', source: 'csvfile', target: 'getSpliters' },
  { id: 'e26', source: 'getSpliters', target: 'target spliter' },
  { id: 'e27', source: 'getSpliters', target: 'pair spliter' },
])
</script>

<template>
  <div class="outer" style="overflow:scroll;">
    <img src="./assets/images/projectLogic.png">
  </div>

  <div style="height: 900px; background: black">
    <VueFlow v-model:nodes="nodes" :edges="edges" :node-types="nodeTypes">
    </VueFlow>
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
/* import the necessary styles for Vue Flow to work */
@import '@vue-flow/core/dist/style.css';

/* import the default theme, this is optional but generally recommended */
@import '@vue-flow/core/dist/theme-default.css';

.vue-flow__node {
    background: green;
    color: white;
    border: 0px solid black;
    border-radius: 4px;
    box-shadow: 0 0 0 4px blue;
    padding: 0px;
}
</style>
