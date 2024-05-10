<script setup>
import { ref } from 'vue';
import axios from 'axios';
const file = ref(null);
const percentage = ref(0);
const props = defineProps({
    'description': String,
    'url': String
});
async function uploadFile(event) {
    file.value = event.target.files;
    try{
        await axios.putForm("/uploadFile/" + props.url, {
            file: file.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <span>{{ description }}</span><br>
        <input type="file" @change="uploadFile">
        <span>{{ percentage }}%</span>
    </div>
</template>